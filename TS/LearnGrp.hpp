#ifndef CMPGRP_HPP
#define CMPGRP_HPP

#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

// The EvoLearn program runs evolutionary simulations
// Copyright (C) 2025  Olof Leimar
// See Readme.md in top repo directory for copyright notice

//**************************** class CompStim *****************************

// This struct represents an instance of a given type of a compound stimulus. An
// instance has a vector x of feature state values x[i] and a vector w_tr of
// contributions to the reward values from the different feature states, such
// that w_tr[i]*x[i] is the expected contribution from stimulus dimension i. The
// actual reward for an agent that selects the compound stimulus is log-normally
// distributed with this expectation. As a special case, it is assumed here that
// a compound stimulus of type cstyp has all x[i] = 0 except x[cstyp] = 1,
// which is the feature that characterises that compound stimulus type

template<typename flt>
struct CompStim {
public:
    using vf_type = std::vector<flt>;    
    CompStim(int a_nsd, flt a_sdR, int a_cstyp = 0) :
        nsd{a_nsd},
        sdR{a_sdR},
        muR{-sdR*sdR/2},
        cstyp{a_cstyp},
        w_tr(nsd, 0),
        x(nsd, 0) { x[cstyp] = 1; }
    flt Reward(flt srn) const { 
        // srn should be a standard normal random number; we assume log-normal
        // variation in rewards
        return std::inner_product(w_tr.begin(), w_tr.end(), x.begin(),
            static_cast<flt>(0)) * std::exp(muR + sdR*srn);
    }
    int cstyp;     // indicates the type of compound stimulus (and its feature)
    int nsd;       // number of stimulus dimensions
    flt sdR;       // log-scale SD of the stochastic variation in reward
    flt muR;       // log-scale mean to compensate for stochastic variation
    vf_type w_tr;  // 'true' expected reward values for the stimulus dimensions
    vf_type x;     // state values (features) for the stimulus dimensions
};


//************************** struct LearnStat *****************************

// This struct stores data on a learning trial, for a learning agent with
// attention

template<typename flt>
struct LearnStat {
public:
    using vf_type = std::vector<flt>;    
    int i_num;    // individual number
    int g_num;    // group number
    int lph;      // learning phase number
    int tstep;    // time step (round or trial) of learning
    int choice;   // choice
    int n_off;    // number of stimulus dimensions with attention turned off
    flt Rew;      // perceived reward
    flt R_ewm;    // ewm reward estimate
    flt hs;       // consistency for chosen stimulus dimension
    flt w_tr;     // true expected value of chosen stimulus dimension
    flt regret;   // regret of chosen stimulus dimension
};


//************************** class LrnGrp *******************************

// This class simulates learning using a meta-learning approach. A learning
// agent experiences a number of bouts or trials. In each bout there is a
// choice between ncs compound stimuli, each characterised by feature states
// x[i] for one or more stimulus dimensions, with i = 1, ..., nsd. The feature
// states denote the absence/presence of a particular feature, indicated as
// 0/1 for the corresponding x[i]. As a special case, it is assumed here that
// a compound stimulus of type cstyp has all x[i] = 0 except x[cstyp] = 1,
// which is the feature that characterises that compound stimulus type (nsd =
// ncs for this special case)

// The time sequence for a trial is, first, that the agent computes the
// estimated reward values of each compound stimulus (based of the estimated
// values of feature states), then the agent selects one of the compound
// stimuli, using the TS procedure, and perceives a reward. Next, the agent
// updates its state (method Update_agent in Phenotyep.hpp), after with there
// is potential mortality from predation.


template<typename PhenType>
class LrnGrp {
public:
    using phen_type = PhenType;
    using vph_type = std::vector<phen_type>;
    using flt = typename phen_type::flt;
    using stat_type = LearnStat<flt>;
    using vs_type = std::vector<stat_type>;
    using comp_stim = CompStim<flt>;
    using vcs_type = std::vector<comp_stim>;
    using vf_type = std::vector<flt>;
    using vvf_type = std::vector<vf_type>;
    using vi_type = std::vector<int>;
    using vvi_type = std::vector<vi_type>;
    using vb_type = std::vector<bool>;
    using vvb_type = std::vector<vb_type>;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_discr = std::discrete_distribution<int>;
    using rand_norm = std::normal_distribution<flt>;
    LrnGrp(int a_g,
        // int a_nsd,
        int a_ncs,
        int a_nph,
        int a_T,
        flt a_mu0,
        flt a_mu1,
        vcs_type& a_vcs,
        const vvf_type& a_fvals,
        const vph_type& a_inds,
        bool a_hist = false) :
        g{a_g},
        // nsd{a_nsd},
        ncs{a_ncs},
        nph{a_nph},
        T{a_T},
        mu0{a_mu0},
        mu1{a_mu1},
        vcs{a_vcs},
        fvals{a_fvals},
        inds{a_inds},
        hist{a_hist}
        { 
            if (hist) {
                stat.reserve(nph*T); 
            }
        }
    vph_type& Get_inds() { return inds; }
    const vs_type& Get_stat() const { return stat; }
    void Learn(rand_eng& eng);

private:
    void Add_stat(const phen_type& agent, int lph, int tstep, 
        flt w_tr = 1, flt regret = 0);
    int g;                 // learning group size
    // int nsd;               // number of stimulus dimensions
    int ncs;               // number of compound stimulus types
    int nph;               // number of phases of learning
    int T;                 // number of learning rounds per phase
    flt mu0;               // mortality parameter
    flt mu1;               // mortality parameter
    vcs_type& vcs;         // compound stimuli
    const vvf_type& fvals; // feature values for compound stimuli
    vph_type inds;         // phenotypes of individuals in the learning group
    bool hist;             // whether to collect learning history
    vs_type stat;          // learning statistics
};

template<typename PhenType>
void LrnGrp<PhenType>::Learn(rand_eng& eng)
{
    // Simulate learning by individuals over nph phases, each of T trials 
    rand_uni uni(0, 1);
    rand_norm nrm(0, 1);    // standard normal random variation
    // alternative 1: run through entire life-time learning, one individual
    // at a time (not suitable for frequency dependence in rewards)
    for (int j = 0; j < g; ++j) {
        phen_type& agent = inds[j];
        // run through learning phases
        for (int l = 0; l < nph; ++l) {
            // set true expected feature values for this phase
            for (int k = 0; k < ncs; ++k) {
                vcs[k].w_tr = fvals[l];
            }
            // record value of best feature for this phase
            flt w_opt = *std::max_element(fvals[l].begin(), fvals[l].end());
            // run through trials in this phase
            for (int trial = 0; trial < T; ++trial) {
                if (agent.alive) {
                    if (agent.n_off < ncs) { // if some attention is on
                        // construct random normals for Thompson sampling
                        vf_type srn(ncs, 0);
                        for (int k = 0; k < ncs; ++k) {
                            srn[k] = nrm(eng);
                        }
                        agent.Make_choice(srn);
                    }
                    int chs = agent.choice;
                    comp_stim& cs = vcs[chs];
                    // get reward from selected compound stimulus
                    agent.Rew = cs.Reward(nrm(eng));
                    // update agent's state
                    agent.Update_agent();
                    // potential mortality 
                    if (uni(eng) < mu0 + mu1 * (ncs - agent.n_off)) {
                        agent.alive = false;
                    }
                    if (hist) {
                        // store learning statistics for this round
                        flt w_tr = cs.w_tr[chs];
                        flt regret = w_opt - w_tr;
                        Add_stat(agent, l, l*T + trial, w_tr, regret);
                    }
                }
            }
        }
    }
}

template<typename PhenType>
void LrnGrp<PhenType>::Add_stat(const phen_type& agent, int lph, int tstep,
    flt w_tr, flt regret)
{
    stat_type st;
    st.i_num = agent.i_num;
    st.g_num = agent.g_num;
    st.lph = lph;
    st.tstep = tstep;
    int chs = agent.choice;
    st.choice = chs;
    st.n_off = agent.n_off;
    st.Rew = agent.Rew;
    st.R_ewm = agent.R_ewm;
    st.hs = agent.hs[chs];
    st.w_tr = w_tr;
    st.regret = regret;
    stat.push_back(st);
}

#endif // CMPGRP_HPP
