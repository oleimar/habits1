#include "toml.hpp"     // to read input parameters from TOML file
#include "rapidcsv.h"   // read gamete data from CSV file
#include "EvoSim.hpp"
#include "hdf5code.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <random>

#ifdef PARA_RUN
#include <omp.h>
#endif

// The EvoLearn program runs evolutionary simulations
// Copyright (C) 2025  Olof Leimar
// See Readme.md in top repo directory for copyright notice

//************************** Get from TOML file **************************

// convenience function to read from TOML input file; this function is used
// to read a value of type T
template<typename T>
void Get(const toml::value& input,
         T& val, const std::string& name)
{
    val = toml::find<T>(input, name);
}


//************************** class EvoInpData ****************************

EvoInpData::EvoInpData(const char* filename) :
      OK(false)
{
    const toml::value idat = toml::parse(filename);

    Get<std::size_t>(idat, max_num_thrds, "max_num_thrds");
    Get<int>(idat, num_loci, "num_loci");
    Get<int>(idat, num_grps, "num_grps");
    Get<int>(idat, G, "G");
    // Get<int>(idat, nsd, "nsd");
    Get<int>(idat, ncs, "ncs");
    Get<int>(idat, nph, "nph");
    Get<int>(idat, T, "T");
    Get<int>(idat, num_yrs, "num_yrs");
    Get<flt>(idat, rho0, "rho0");
    Get<flt>(idat, sd0, "sd0");
    Get<flt>(idat, sdR, "sdR");
    Get<flt>(idat, alph_hs, "alph_hs");
    Get<flt>(idat, mu0, "mu0");
    Get<flt>(idat, mu1, "mu1");
    Get<std::string>(idat, CSName, "CSName");

    Get<bool>(idat, samp_fvals, "samp_fvals");
    Get<bool>(idat, hist, "hist");
    Get<bool>(idat, read_from_file, "read_from_file");
    Get<std::string>(idat, gam_data_csv, "gam_data_csv");
    if (read_from_file) {
        Get<std::string>(idat, h5InName, "h5InName");
    }
    Get<std::string>(idat, h5OutName, "h5OutName");
    if (hist) {
        Get<std::string>(idat, h5HistName, "h5HistName");
    }

    InpName = std::string(filename);
    OK = true;
}


//****************************** Class Evo *****************************

Evo::Evo(const EvoInpData& eid) :
    id{eid},
    num_loci{id.num_loci},
    num_grps{id.num_grps},
    G{id.G},
    N{num_grps*G},
    // nsd{id.nsd},
    ncs{id.ncs},
    nph{id.nph},
    T{id.T},
    num_yrs{id.num_yrs},
    rho0{id.rho0},
    sd0{id.sd0},
    sdR{id.sdR},
    alph_hs{id.alph_hs},
    mu0{id.mu0},
    mu1{id.mu1},
    hist{id.hist},
    num_thrds{1},
    fvals(nph, std::vector<flt>(ncs, 0))
{
    // decide on number of threads for parallel processing
#ifdef PARA_RUN
    num_thrds = omp_get_max_threads();
    if (num_thrds > id.max_num_thrds) num_thrds = id.max_num_thrds;
    // check that there is at least one local group per thread
    if (num_thrds > num_grps) num_thrds = num_grps;
    std::cout << "Number of threads: "
              << num_thrds << '\n';
#endif
    // generate one seed for each thread
    sds.resize(num_thrds);
    std::random_device rd;
    for (unsigned i = 0; i < num_thrds; ++i) {
        // set up thread-local to be random number engines
        sds[i] = rd();
        vre.push_back(rand_eng(sds[i]));
    }
    // read gamete data from CSV file
    rapidcsv::Document gam_data(id.gam_data_csv);
    // Probability of mutation at each locus
    vf_type mut_rate = gam_data.GetColumn<flt>(1);
    // SD of mutational increments at each locus
    vf_type SD = gam_data.GetColumn<flt>(2);
    // Maximal allelic value at each locus
    vf_type max_val = gam_data.GetColumn<flt>(3);
    // Minimal allelic value at each locus
    vf_type min_val = gam_data.GetColumn<flt>(4);
    // Recombination rates
    vf_type rho = gam_data.GetColumn<flt>(5);
    for (unsigned i = 0; i < num_thrds; ++i) {
        // set up thread-local to be mutation records, with thread-local engine
        // and parameters controlling mutation, segregation and recombination
        rand_eng& eng = vre[i];
        mut_rec_type mr(eng, num_loci);
        for (unsigned l = 0; l < num_loci; ++l) {
            mr.mut_rate[l] = mut_rate[l];
            mr.SD[l] = SD[l];
            mr.max_val[l] = max_val[l];
            mr.min_val[l] = min_val[l];
            mr.rho[l] = rho[l];
        }
        vmr.push_back(mr);
    }

    int thread_num = 0;

    // construct global container of compound stimuli
    comp_stim cs0(ncs, sdR);
    vcs.resize(ncs, cs0);
    // read data for compound stimuli from CSV file
    rapidcsv::Document CSdat(id.CSName);
    for (int k = 0; k < ncs; ++k) {
        comp_stim& cs = vcs[k];
        cs.cstyp = k;
        cs.nsd = ncs;
        cs.sdR = sdR;
        // CSdat contains features for the compound stimuli in the 
        // first ncs columns: x[k] = 1 for cstyp = k
        cs.x = CSdat.GetColumn<flt>(k);
    }
    if (!id.samp_fvals) {
        // extract values of stimulus features in different phases from CSV
        // file (this becomes the same as values for the compound stimuli)
        for (int l = 0; l < nph; ++l) {
            fvals[l] = CSdat.GetColumn<flt>(ncs + l);
        }
    }

    // Note concerning thread safety: in order to avoid possible problems with
    // multiple threads, the std::vector containers pop, and stat are
    // allocated once and for all here, and thread-local data are then copied
    // into position in these (thus avoiding potentially unsafe push_back and
    // insert).

    // Create N "placeholder individuals" in the population (based on the
    // constructor these are not alive). The convention is that individuals in
    // local group grp (grp = 0, ..., num_grps-1) are found as alive
    // phenotypes at pop[j] with j = grp*G, Grp*G + 1, ..., grp*G + G - 1,
    // i.e. j = grp*G + i with i = 0, ..., G-1. These values of i (i.e., i =
    // 0, ..., G) correspond to the original values of i_num in each group
    // g_num = grp, assigned to individuals when read_from_file is false.
    gam_type gam(num_loci);
    ind_type ind(ncs, gam);
    // set phenotype parameters
    phen_type& ph = ind.phenotype;
    ph.alph_hs = alph_hs;
    pop.resize(N, ind);
    // create offspring placeholder genotypes
    gen_type gen(gam);
    offspr.resize(N, gen);
    // history stats
    if (hist) {
        stat.reserve(N*T*nph);
    }

    // check if population data should be read from file
    if (id.read_from_file) {
        // Read_pop(id.InName);
        h5_read_pop(id.h5InName);
        // form new population from offspring
        ReplacePop();
    } else {
        // construct all individuals with the same genotypes
        gam_type gam(num_loci); // starting gamete
        // Starting allelic values
        vf_type all0 = gam_data.GetColumn<flt>(6);
        for (int l = 0; l < num_loci; ++l) {
            gam.gam_dat[l] = all0[l];
        }
        for (int grp = 0; grp < num_grps; ++grp) { // local groups
            for (int i = 0; i < G; ++i) { // individuals in group
                int j = grp*G + i;
                ind_type ind(ncs, gam);
                phen_type& ph = ind.phenotype;
                ph.alph_hs = alph_hs;
                ph.i_num = i;      // set individual number in group
                ph.g_num = grp;    // set local group number
                ph.female = true;
                ph.alive = true;
                pop[j] = ind;
            }
        }
    }
}


void Evo::Run()
{
    Timer timer(std::cout);
    timer.Start();
    ProgressBar PrBar(std::cout, num_yrs);
    // Run through years. Time sequence within a year:

    // (1) set up local groups; (2) learning in each group, including
    // mortality during foraging; (3) reproduction, with parents of a new
    // individual selected fro entire population, with probabilities
    // proportional to life-time accumulated rewards; (4) form population for
    // next year from the offspring.
    for (unsigned yr = 0; yr < num_yrs; ++yr) {
        // use parallel for processing over the learning in local
        // groups, within a year
#pragma omp parallel for num_threads(num_thrds)
        for (unsigned grp = 0; grp < num_grps; ++grp) {
#ifdef PARA_RUN
            int thread_num = omp_get_thread_num();
#else
            int thread_num = 0;
#endif
            // thread-local random number engine
            rand_eng& eng = vre[thread_num];
            // thread-local mutation record
            mut_rec_type& mr = vmr[thread_num];
            // distributions needed
            rand_uni uni(0, 1);
            // (1) set up local groups
            // thread-local container for local group phenotypes
            vph_type loc_ph;
            loc_ph.reserve(G);
            for (int i = 0; i < G; ++i) {
                int j = grp*G + i;
                phen_type& ph = pop[j].phenotype;
                loc_ph.push_back(ph);
            }
            // thread-local container for compound stimuli; needs to be local
            // because true expected values of features can change during
            // learning in a group
            vcs_type loc_vcs = vcs;
            // thread-local container for feature values; needs to be local
            // because true expected values can be randomly drawn
            vvf_type loc_fvals = fvals;
            if (id.samp_fvals) {
                // sample the feature values in the different phases from
                // a log-scale normal prior with mean rho0 and SD of sd0
                rand_norm prio_log_scale(rho0, sd0); 
                for (int l = 0; l < nph; ++l) {
                    for (int k = 0; k < ncs; ++k) {
                        flt rho = prio_log_scale(eng);
                        loc_fvals[l][k] = std::exp(rho);
                    }
                }
            }
            // get history only for single threaded and final year
            bool hist1 = false;
            if (num_thrds == 1 && yr == num_yrs - 1 && hist) {
                hist1 = true;
            }
            // (2) learning in each group
            lg_type lg(G, ncs, nph, T, mu0, mu1, loc_vcs,
                loc_fvals, loc_ph, hist1);
            lg.Learn(eng);
            vph_type& lph = lg.Get_inds();
            for (int i = 0; i < G; ++i) {
                loc_ph[i] = lph[i];
            }
            // copy phenotypes from loc_ph into population
            for (unsigned i = 0; i < G; ++i) {
                int j = grp*G + i;
                pop[j].phenotype = loc_ph[i];
            }
            if (hist1) {
                // append history
                const vs_type& st = lg.Get_stat();
                stat.insert(stat.end(), st.begin(), st.end());
            }
        } // end of parallel for (over local groups)

        // (3) reproduction
        // reproduction in entire population
        SelectReproduce(vmr[0]);
        // if not final year, get new pop from offspring
        if (yr < num_yrs - 1) {
            ReplacePop();
        }
        // all set to start next year
        ++PrBar;
    }

    PrBar.Final();
    timer.Stop();
    timer.Display();
    h5_write_pop(id.h5OutName);
    if (hist) {
        h5_write_hist(id.h5HistName);
    }
}

void Evo::SelectReproduce(mut_rec_type& mr)
{
    // get discrete distribution with lifetime rewards as weights
    vf_type wei(N);
    flt tot_wei = 0;
    for (int grp = 0; grp < num_grps; ++grp) {
        for (int i = 0; i < G; ++i) {
            int j = grp*G + i;
            const phen_type& ph = pop[j].phenotype;
            wei[j] = ph.R_tot;
            tot_wei += wei[j];
        }
    }
    if (tot_wei > 0) {
        // assume that individuals are hermaphrodites
        // get discrete distribution for parents
        rand_discr dscr(wei.begin(), wei.end());
        // get offspring genotypes
        for (int grp = 0; grp < num_grps; ++grp) {
            for (int i = 0; i < G; ++i) {
                int j = grp*G + i;
                // find "mother" for individual to be constructed
                int imat = dscr(mr.eng);
                const ind_type& matind = pop[imat];
                // find "father" for individual to be constructed
                int ipat = dscr(mr.eng);
                const ind_type& patind = pop[ipat];
                // form new genotype and copy to offspr
                gen_type gen(matind.GetGamete(mr), patind.GetGamete(mr));
                offspr[j] = gen;
            }
        }
    } else { // no reproduction (should not happen)
        std::cout << "NOTE: reproduction failed\n";
    }
}

void Evo::ReplacePop()
{
    for (int grp = 0; grp < num_grps; ++grp) { // local groups
        for (int i = 0; i < G; ++i) { // individuals in group
            int j = grp*G + i;
            // replace individuals with offspring from genotypes
            gen_type& gen = offspr[j];
            pop[j].genotype = gen;
            phen_type& ph = pop[j].phenotype;
            ph.Assign(gen, alph_hs);
            // assign data not set correctly by Assign
            ph.i_num = i;      // set individual number in group
            ph.g_num = grp;    // set local group number
            ph.female = true;
            ph.alive = true;
        }
    }
}

// Helper function for reading and assigning gamete data
void Evo::read_and_assign_gam(h5R& h5, const std::string& dataset_name,
    vvf_type& gams, vgen_type& vgen, gam_type gen_type::*member) {
    h5.read_flt_arr(dataset_name, gams);
    for (int i = 0; i < gams.size(); ++i) {
        gen_type& gen = vgen[i];
        gam_type& gam = gen.*member;
        for (int l = 0; l < num_loci; ++l) {
            gam[l] = gams[i][l];
        }
    }
}

// Helper function to read and assign flt attributes
void Evo::read_and_assign_flt(h5R& h5, const std::string& dataset_name,
    vf_type& fval, flt phen_type::*member) {
    h5.read_flt(dataset_name, fval);
    for (int i = 0; i < fval.size(); ++i) {
        pop[i].phenotype.*member = fval[i];
    }
}

// Helper function to read and assign int attributes
void Evo::read_and_assign_int(h5R& h5, const std::string& dataset_name,
    vi_type& i_val, int phen_type::*member) {
    h5.read_int(dataset_name, i_val);
    for (int i = 0; i < i_val.size(); ++i) {
        pop[i].phenotype.*member = i_val[i];
    }
}

// Helper function to read and assign bool attributes
void Evo::read_and_assign_bool(h5R& h5, const std::string& dataset_name,
    std::vector<int>& ival, bool phen_type::*member) {
    h5.read_int(dataset_name, ival);
    for (int i = 0; i < ival.size(); ++i) {
        pop[i].phenotype.*member = ival[i];
    }
}

// Helper function for reading and assigning vf_type data
void Evo::read_and_assign_vf_type(h5R& h5, const std::string& dataset_name,
    vvf_type& flt_pars, vf_type phen_type::*member) {
    h5.read_flt_arr(dataset_name, flt_pars);
    for (int i = 0; i < flt_pars.size(); ++i) {
        phen_type& ph = pop[i].phenotype;
        for (int j = 0; j < flt_pars[i].size(); ++j) {
            (ph.*member)[j] = flt_pars[i][j];
        }
    }
}

// Helper function for reading and assigning vi_type data
void Evo::read_and_assign_vi_type(h5R& h5, const std::string& dataset_name,
    vvi_type& i_pars, vi_type phen_type::*member) {
    h5.read_int_arr(dataset_name, i_pars);
    for (int i = 0; i < i_pars.size(); ++i) {
        phen_type& ph = pop[i].phenotype;
        for (int j = 0; j < i_pars[i].size(); ++j) {
            (ph.*member)[j] = i_pars[i][j];
        }
    }
}

// Helper function to write phenotype flt data
void Evo::write_flt_attr(h5W& h5, const std::string& dataset_name,
    vf_type& attr_values, flt phen_type::*member) const {
    for (int i = 0; i < attr_values.size(); ++i) {
        attr_values[i] = pop[i].phenotype.*member;
    }
    h5.write_flt(dataset_name, attr_values);
}

// Helper function to write phenotype int data
void Evo::write_int_attr(h5W& h5, const std::string& dataset_name,
    vi_type& attr_values, int phen_type::*member) const {
    for (int i = 0; i < attr_values.size(); ++i) {
        attr_values[i] = pop[i].phenotype.*member;
    }
    h5.write_int(dataset_name, attr_values);
}

// Helper function to write phenotype bool data
void Evo::write_bool_attr(h5W& h5, const std::string& dataset_name,
    std::vector<int>& attr_values, bool phen_type::*member) const {
    for (int i = 0; i < attr_values.size(); ++i) {
        attr_values[i] = pop[i].phenotype.*member;
    }
    h5.write_int(dataset_name, attr_values);
}

// Helper function to write vf_type data
void Evo::write_vf_type(h5W& h5, const std::string& dataset_name,
    vvf_type& flt_pars, vf_type phen_type::*member) const {
    for (int i = 0; i < flt_pars.size(); ++i) {
        const phen_type& ph = pop[i].phenotype;
        for (int j = 0; j < flt_pars[i].size(); ++j) {
            flt_pars[i][j] = (ph.*member)[j];
        }
    }
    h5.write_flt_arr(dataset_name, flt_pars);
}

// Helper function to write vi_type data
void Evo::write_vi_type(h5W& h5, const std::string& dataset_name,
    vvi_type& i_pars, vi_type phen_type::*member) const {
    for (int i = 0; i < i_pars.size(); ++i) {
        const phen_type& ph = pop[i].phenotype;
        for (int j = 0; j < i_pars[i].size(); ++j) {
            i_pars[i][j] = (ph.*member)[j];
        }
    }
    h5.write_int_arr(dataset_name, i_pars);
}

// Helper function to write stat_type flt data
void Evo::write_stat_flt(h5W& h5, const std::string& dataset_name,
    vf_type& attr_values, 
    flt stat_type::*member) const {
    for (int i = 0; i < attr_values.size(); ++i) {
        attr_values[i] = stat[i].*member;
    }
    h5.write_flt(dataset_name, attr_values);
}

// Helper function to write stat_type int data
void Evo::write_stat_int(h5W& h5, const std::string& dataset_name,
    vi_type& attr_values, int stat_type::*member) const {
    for (int i = 0; i < attr_values.size(); ++i) {
        attr_values[i] = stat[i].*member;
    }
    h5.write_int(dataset_name, attr_values);
}

// Helper function to write stat_type vf_type data
void Evo::write_stat_vf_type(h5W& h5, const std::string& dataset_name,
    vvf_type& attr_values, vf_type stat_type::*member) const {
    for (int i = 0; i < attr_values.size(); ++i) {
        const stat_type& st = stat[i];
        for (int j = 0; j < attr_values[i].size(); ++j) {
            attr_values[i][j] = (st.*member)[j];
        }
    }
    h5.write_flt_arr(dataset_name, attr_values);
}

void Evo::h5_read_pop(const std::string& infilename)
{
    // read data and put in pop and offspr
    h5R h5(infilename);
    // gametes of individuals in pop
    vvf_type gams(N, vf_type(num_loci));
    h5.read_flt_arr("MatGam", gams);
    for (int i = 0; i < N; ++i) {
        gam_type& gam = pop[i].genotype.mat_gam;
        for (int l = 0; l < num_loci; ++l) {
            gam[l] = gams[i][l];
        }
    }
    h5.read_flt_arr("PatGam", gams);
    for (int i = 0; i < N; ++i) {
        gam_type& gam = pop[i].genotype.pat_gam;
        for (int l = 0; l < num_loci; ++l) {
            gam[l] = gams[i][l];
        }
    }
    // Read offspring gametes
    read_and_assign_gam(h5, "OffsMatGam", gams, offspr,
        &gen_type::mat_gam);
    read_and_assign_gam(h5, "OffsPatGam", gams, offspr,
        &gen_type::pat_gam);
    // Read phenotype attributes
    vf_type f_val(N);
    // flt attributes
    read_and_assign_flt(h5, "w0", f_val, &phen_type::w0);
    read_and_assign_flt(h5, "alph_w", f_val, &phen_type::alph_w);
    read_and_assign_flt(h5, "lambda", f_val, &phen_type::lambda);
    read_and_assign_flt(h5, "thr1", f_val, &phen_type::thr1);
    read_and_assign_flt(h5, "thr2", f_val, &phen_type::thr2);
    read_and_assign_flt(h5, "t_thr", f_val, &phen_type::t_thr);
    read_and_assign_flt(h5, "beta", f_val, &phen_type::beta);
    read_and_assign_flt(h5, "Rew", f_val, &phen_type::Rew);
    read_and_assign_flt(h5, "R_tot", f_val, &phen_type::R_tot);
    read_and_assign_flt(h5, "delt", f_val, &phen_type::delt);
    // int attributes
    vi_type i_val(N);
    read_and_assign_int(h5, "choice", i_val, &phen_type::choice);
    read_and_assign_int(h5, "n_Rew", i_val, &phen_type::n_Rew);
    read_and_assign_int(h5, "n_off", i_val, &phen_type::n_off);
    read_and_assign_int(h5, "i_num", i_val, &phen_type::i_num);
    read_and_assign_int(h5, "g_num", i_val, &phen_type::g_num);
    // bool attributes
    std::vector<int> b_val(N);
    read_and_assign_bool(h5, "female", b_val, &phen_type::female);
    read_and_assign_bool(h5, "alive", b_val, &phen_type::alive);
    // vf_type attributes
    vvf_type flt_pars(N, vf_type(ncs));
    read_and_assign_vf_type(h5, "w", flt_pars, &phen_type::w);
    read_and_assign_vf_type(h5, "att", flt_pars, &phen_type::att);
    vvi_type int_pars(N, vi_type(ncs));
    read_and_assign_vi_type(h5, "n_feat", int_pars, &phen_type::n_feat);
    read_and_assign_vi_type(h5, "l_feat", int_pars, &phen_type::l_feat);
}

void Evo::h5_write_pop(const std::string& outfilename) const
{
    h5W h5(outfilename);
    // gametes
    std::vector<vf_type> gams(N, vf_type(num_loci));
    // write maternal gametes
    for (int i = 0; i < N; ++i) {
        const gam_type& gam = pop[i].genotype.mat_gam;
        for (int l = 0; l < num_loci; ++l) {
            gams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("MatGam",gams);
    // write paternal gametes
    for (int i = 0; i < N; ++i) {
        const gam_type& gam = pop[i].genotype.pat_gam;
        for (int l = 0; l < num_loci; ++l) {
        gams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("PatGam",gams);
    // offspring gametes
    // write maternal gametes of offspring
    for (unsigned i = 0; i < N; ++i) {
        const gam_type& gam = offspr[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("OffsMatGam", gams);
    // write paternal gametes of offspring
    for (unsigned i = 0; i < N; ++i) {
        const gam_type& gam = offspr[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("OffsPatGam", gams);

    // write flt phenotype attributes
    vf_type f_val(N);
    write_flt_attr(h5, "w0", f_val, &phen_type::w0);
    write_flt_attr(h5, "alph_w", f_val, &phen_type::alph_w);
    write_flt_attr(h5, "lambda", f_val, &phen_type::lambda);
    write_flt_attr(h5, "thr1", f_val, &phen_type::thr1);
    write_flt_attr(h5, "thr2", f_val, &phen_type::thr2);
    write_flt_attr(h5, "t_thr", f_val, &phen_type::t_thr);
    write_flt_attr(h5, "beta", f_val, &phen_type::beta);
    write_flt_attr(h5, "Rew", f_val, &phen_type::Rew);
    write_flt_attr(h5, "R_ewm", f_val, &phen_type::R_ewm);
    write_flt_attr(h5, "R_tot", f_val, &phen_type::R_tot);
    write_flt_attr(h5, "delt", f_val, &phen_type::delt);
    // write int phenotype attributes
    vi_type i_val(N);
    write_int_attr(h5, "choice", i_val, &phen_type::choice);
    write_int_attr(h5, "n_Rew", i_val, &phen_type::n_Rew);
    write_int_attr(h5, "n_off", i_val, &phen_type::n_off);
    write_int_attr(h5, "i_num", i_val, &phen_type::i_num);
    write_int_attr(h5, "g_num", i_val, &phen_type::g_num);
    // write bool phenotype attributes
    std::vector<int> b_val(N);
    write_bool_attr(h5, "female", b_val, &phen_type::female);
    write_bool_attr(h5, "alive", b_val, &phen_type::alive);
    // write vf_type phenotype attributes
    vvf_type flt_pars(N, vf_type(ncs));
    write_vf_type(h5, "w", flt_pars, &phen_type::w);
    write_vf_type(h5, "att", flt_pars, &phen_type::att);
    write_vf_type(h5, "hs", flt_pars, &phen_type::hs);
    // write vi_type phenotype attributes
    vvi_type int_pars(N, vi_type(ncs));
    write_vi_type(h5, "n_feat", int_pars, &phen_type::n_feat);
    write_vi_type(h5, "l_feat", int_pars, &phen_type::l_feat);
}

void Evo::h5_write_hist(const std::string& histfilename) const
{
    h5W h5(histfilename);
    unsigned hlen = stat.size();
    // std::vector to hold int data
    vi_type i_val(hlen);
    write_stat_int(h5, "i_num", i_val, &stat_type::i_num);
    write_stat_int(h5, "g_num", i_val, &stat_type::g_num);
    write_stat_int(h5, "lph", i_val, &stat_type::lph);
    write_stat_int(h5, "tstep", i_val, &stat_type::tstep);
    write_stat_int(h5, "choice", i_val, &stat_type::choice);
    write_stat_int(h5, "n_off", i_val, &stat_type::n_off);
    // std::vector to hold flt data
    vf_type f_val(hlen);
    write_stat_flt(h5, "Rew", f_val, &stat_type::Rew);
    write_stat_flt(h5, "delt", f_val, &stat_type::delt);
    write_stat_flt(h5, "R_ewm", f_val, &stat_type::R_ewm);
    write_stat_flt(h5, "w", f_val, &stat_type::w);
    write_stat_flt(h5, "hs", f_val, &stat_type::hs);
    write_stat_flt(h5, "w_tr", f_val, &stat_type::w_tr);
    write_stat_flt(h5, "regret", f_val, &stat_type::regret);
}
