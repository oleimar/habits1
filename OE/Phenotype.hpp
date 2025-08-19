#ifndef PHENOTYPE_HPP
#define PHENOTYPE_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

// The EvoLearn program runs evolutionary simulations
// Copyright (C) 2025  Olof Leimar
// See Readme.md in top repo directory for copyright notice

//************************* struct Phenotype ******************************

template<typename GenType>
struct Phenotype {
// public:
    using flt = float;
    using vf_type = std::vector<flt>;
    using vi_type = std::vector<int>;
    using gen_type = GenType;
    using val_type = typename gen_type::val_type;
    Phenotype(int a_ncs = 2,
              const gen_type& gt = gen_type()) :
        ncs{a_ncs},
        w(ncs, 0),
        att(ncs, 1),
        hs(ncs, 0),
        n_feat(ncs, 0),
        l_feat(ncs, 0)
        { Assign(gt); }
    void Assign(const gen_type& gt, 
        flt a_alph_hs = 0.25);
    // Update of agent state
    void Update_agent();
    bool Female() const { return female; }
    // public data members
    flt w0;       // starting estimate of reward value
    flt alph_w;   // Rescorla-Wagner learning rate for w
    flt lambda;   // ewm factor for estimate of recent rewards
    flt thr1;     // threshold for turning off attention
    flt thr2;     // threshold for turning on attention
    flt t_thr;    // threshold for number of rounds since attention was off
    flt beta;     // soft-max choice parameter (inverse temperature)
    flt alph_hs;  // habit strength learning rate
    flt Rew;      // perceived reward
    flt R_ewm;    // estimate of recent rewards
    flt R_tot;    // lifetime reward (so far)
    flt delt;     // prediction error (TD error)
    int ncs;      // number of stimulus dimensions (= number of compound stims)
    int choice;   // choice
    int n_Rew;    // number of rewards so far since attention was off
    int n_off;    // number of stimulus dimensions with attention off
    int i_num;    // individual number in local group
    int g_num;    // local group number
    bool female;
    bool alive;
    vf_type w;         // estimates of reward values
    vf_type att;       // attention (1/0) to stimulus dimensions
    vf_type hs;        // consistency of choice (habit strength) estimate
    vi_type n_feat;    // number of times feature is observed on chosen CS
    vi_type l_feat;    // number of rounds lacking feature on chosen CS
};

template<typename GenType>
void Phenotype<GenType>::Assign(const gen_type& gt, 
    flt a_alph_hs)
{
    val_type val = gt.Value();
    // assume val is a vector with components corresponding to the traits 
    // w0, alph_w, lambda, thr1, thr2, t_thr, beta
    w0 = val[0];
    alph_w = val[1];
    lambda = val[2];
    thr1 = val[3];
    thr2 = val[4];
    t_thr = val[5];
    beta = val[6];
    alph_hs = a_alph_hs;
    Rew = 0;
    R_ewm = 0;
    R_tot = 0;
    delt = 0;
    choice = 0;
    n_Rew = 0;
    n_off = 0;
    i_num = 0;
    g_num = 0;
    female = false;
    alive = false;
    for (int k = 0; k < ncs; ++k) {
        w[k] = w0;
        att[k] = 1;
        hs[k] = 0;
        n_feat[k] = 0;
        l_feat[k] = 0;
    }
}

// Meta-learning update of agent state. The agent first updates the
// consistency of choice (for the chosen option), then updates the ewm
// estimate of rewards, then decides whether to turn off attention. The agent
// then checks if recent rewards have become so small that all attention
// should be turned on. Finally, for the dimensions attended to, the agent
// updates the estimated reward value w[i] for the feature of the selected
// compound stimulus.

template<typename GenType>
void Phenotype<GenType>::Update_agent() {
    // update consistency of choice
    for (int k = 0; k < ncs; ++k) {
        if (att[k] > 0) {
            if (choice == k) {
                hs[k] += alph_hs*(1 - hs[k]);
            } else {
                hs[k] += alph_hs*(-hs[k]);
            }
        }
    }
    // update R_ewm estimate of recent rewards
    R_ewm += lambda*(Rew - R_ewm);
    ++n_Rew;
    // prediction error
    delt = Rew - w[choice];
    // update book-keeping counts
    for (int k = 0; k < ncs; ++k) {
        if (att[k] > 0) {
            if (choice == k) {
                l_feat[k] = 0;
                ++n_feat[k];
            } else {
                // count rounds since feature was last chosen
                ++l_feat[k];
            }
        }
    }
    // check if attention to dimensions should be turned off
    for (int k = 0; k < ncs; ++k) {
        if (att[k] > 0) {
            if (choice == k) {
                if (hs[k] > thr1) {
                    // this is a "trigger" to turn off all attention
                    for (int j = 0; j < ncs; ++j) {
                        if (att[j] > 0) {
                            att[j] = 0;
                            n_feat[j] = 0;
                            l_feat[j] = 0;
                            ++n_off;
                        }
                    }
                }
            } else if (l_feat[k] > t_thr) {
                att[k] = 0;
                n_feat[k] = 0;
                l_feat[k] = 0;
                ++n_off;
            }
        }
    }
    // check if all attention should be turned on (with reset of state)
    if (n_off > 0 && n_Rew > t_thr && R_ewm < thr2) {
        // turn on attention and reset for all dimensions
        for (int k = 0; k < ncs; ++k) {
            att[k] = 1;
            // reset estimated reward values
            w[k] = w0;
            // reset consistency of choice to zero
            hs[k] = 0;
            // reset count of feature observations
            n_feat[k] = 0;
            l_feat[k] = 0;
        }
        // attention is on for all dimensions
        n_off = 0;
        // reset count of rewards
        n_Rew = 0;
    }
    // update estimated value of chosen feature
    if (att[choice] > 0) w[choice] += alph_w*delt;
    // update total rewards (expected reproductive success)
    R_tot += Rew;
}

#endif // PHENOTYPE_HPP
