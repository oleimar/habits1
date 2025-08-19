#ifndef EVOSIM_HPP
#define EVOSIM_HPP

// NOTE: for easier debugging one can comment out the following
#ifdef _OPENMP
#define PARA_RUN
#endif

#include "Genotype.hpp"
#include "Phenotype.hpp"
#include "Individual.hpp"
#include "LearnGrp.hpp"
#include "hdf5code.hpp"
#include <vector>
#include <string>
#include <random>
#include <cmath>

// The EvoLearn program runs evolutionary simulations
// Copyright (C) 2025  Olof Leimar
// See Readme.md in top repo directory for copyright notice

// An individual has 7 genetically determined traits w0, alph_w, lambda, thr1,
// thr2, t_thr, beta (see Phenotype.hpp) and there is one locus for each trait

//************************* Class EvoInpData ***************************

// This class is used to 'package' input data in a single place; the
// constructor extracts data from an input file

class EvoInpData {
public:
    // the type flt should be the same type as defined in struct Phenotype
    using flt =float;
    using vf_type = std::vector<flt>;
    // using vint_type = std::vector<int_type>;
    std::size_t max_num_thrds;  // Max number of threads to use
    int num_loci;               // Number of loci in individual's genotype
    int num_grps;               // Number of local groups
    int G;                      // Local group size
    // int nsd;                    // number of stimulus dimensions (=ncs)
    int ncs;                    // number of compound stimuli
    int nph;                    // number of phases of learning
    int T;                      // number of learning rounds each phase
    int num_yrs;                // number of years to simulate
    flt rho0;                   // mean of log-scale prior feature values
    flt sd0;                    // sd of log-scale prior feature values
    flt sdR;                    // sd of log-scale reward variation
    flt alph_hs;                // habit strength learning rate
    flt mu0;                    // mortality parameter
    flt mu1;                    // mortality parameter
    std::string CSName;         // file name for CSV file for comp stims
    bool samp_fvals;            // Whether to randomly sample feature values
    bool hist;                  // Whether to compute and save history
    bool read_from_file;        // Whether to read population from file
    std::string gam_data_csv;   // File name for input of gamete data
    std::string h5InName;       // File name for input of population
    std::string h5OutName;      // File name for output of population
    std::string h5HistName;     // File name for output of history

    std::string InpName;  // Name of input data file
    bool OK;              // Whether input data has been successfully read

    EvoInpData(const char* filename);
};


//***************************** Class Evo ******************************

class Evo {
public:
    // types needed to define individual
    using mut_rec_type = MutRec<MutIncrNorm<>>;
    using gam_type = Gamete<mut_rec_type>;
    using gen_type = Diplotype<gam_type>;
    using phen_type = Phenotype<gen_type>;
    using ind_type = Individual<gen_type, phen_type>;
    using flt = typename phen_type::flt;
    using stat_type = LearnStat<flt>;
    using comp_stim = CompStim<flt>;
    // use std::vector container for population
    using vind_type = std::vector<ind_type>;
    using vgen_type = std::vector<gen_type>;
    using vph_type = std::vector<phen_type>;
    using lg_type = LrnGrp<phen_type>;
    using vf_type = std::vector<flt>;
    using vvf_type = std::vector<vf_type>;
    using vi_type = std::vector<int>;
    using vvi_type = std::vector<vi_type>;
    using vcs_type = std::vector<comp_stim>;
    using vs_type = std::vector<stat_type>;
    using vui_type = std::vector<unsigned>;
    using rand_eng = std::mt19937;
    using vre_type = std::vector<rand_eng>;
    using vmr_type = std::vector<mut_rec_type>;
    using rand_u = std::uniform_int_distribution<int>;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_norm = std::normal_distribution<flt>;
    using rand_discr = std::discrete_distribution<int>;
    Evo(const EvoInpData& eid);
    void Run();
    void h5_read_pop(const std::string& infilename);
    void h5_write_pop(const std::string& outfilename) const;
    void h5_write_hist(const std::string& histfilename) const;
private:
    void SelectReproduce(mut_rec_type& mr);
    void ReplacePop();
    // Helper function for reading and assigning gamete data
    void read_and_assign_gam(h5R& h5, const std::string& dataset_name,
        vvf_type& gams, vgen_type& vgen, gam_type gen_type::*member);
    // Helper functions for reading and assigning phenotype data
    void read_and_assign_flt(h5R& h5, const std::string& dataset_name,
        vf_type& fval, flt phen_type::*member);
    void read_and_assign_int(h5R& h5, const std::string& dataset_name,
        vi_type& i_val, int phen_type::*member);
    void read_and_assign_bool(h5R& h5, const std::string& dataset_name,
        std::vector<int>& ival, bool phen_type::*member);
    // Helper function for reading and assigning vf_type data
    void read_and_assign_vf_type(h5R& h5, const std::string& dataset_name,
        vvf_type& flt_pars, vf_type phen_type::*member);
    // Helper function for reading and assigning vi_type data
    void read_and_assign_vi_type(h5R& h5, const std::string& dataset_name,
        vvi_type& int_pars, vi_type phen_type::*member);
    // Helper function to write phenotype flt data
    void write_flt_attr(h5W& h5, const std::string& dataset_name,
                        vf_type& attr_values, flt phen_type::*member) const;
    // Helper function to write phenotype int data
    void write_int_attr(h5W& h5, const std::string& dataset_name,
                         vi_type& attr_values, 
                         int phen_type::*member) const;
    // Helper function to write phenotype bool data
    void write_bool_attr(h5W& h5, const std::string& dataset_name,
                         std::vector<int>& attr_values, 
                         bool phen_type::*member) const;
    // Helper function to write vf_type data
    void write_vf_type(h5W& h5, const std::string& dataset_name,
                       vvf_type& attr_values, 
                       vf_type phen_type::*member) const;
    // Helper function to write vi_type data
    void write_vi_type(h5W& h5, const std::string& dataset_name,
                        vvi_type& attr_values, 
                        vi_type phen_type::*member) const;
    // Helper function to write stat_type flt data
    void write_stat_flt(h5W& h5, const std::string& dataset_name,
                        vf_type& attr_values, flt stat_type::*member) const;
    // Helper function to write stat_type int data
    void write_stat_int(h5W& h5, const std::string& dataset_name,
                        vi_type& attr_values, int stat_type::*member) const;
    // Helper function to write stat_type vf_type data
    void write_stat_vf_type(h5W& h5, const std::string& dataset_name,
                        vvf_type& attr_values, vf_type stat_type::*member) const;


    EvoInpData id;
    int num_loci;               // Number of loci in individual's genotype
    int num_grps;               // Number of local groups
    int G;                      // Local group size
    int N;                      // Population size
    // int nsd;                    // number of stimulus dimensions
    int ncs;                    // number of compound stimulus types to use
    int nph;                    // number of phases of learning
    int T;                      // number of learning trials per replicate
    int num_yrs;                // number of years to simulate
    flt rho0;                   // mean of log-scale prior feature values
    flt sd0;                    // sd of log-scale prior feature values
    flt sdR;                    // sd of log-scale reward variation
    flt alph_hs;                // habit strength learning rate
    flt mu0;                    // mortality parameter
    flt mu1;                    // mortality parameter
    bool hist;
    std::size_t num_thrds;
    vcs_type vcs;               // vector of compound stimuli
    vvf_type fvals;             // vector of feature values
    vui_type sds;               // seeds for random numbers
    vre_type vre;               // random number engines
    vmr_type vmr;               // mutation records
    vind_type pop;              // population container
    vgen_type offspr;           // container for offspring genotypes
    vs_type stat;               // learning statistics
};

#endif // EVOSIM_HPP
