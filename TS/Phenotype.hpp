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

// Assumptions about GenType:
// types:
//   val_type
// member functions:
//   val_type Value()

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
        y_bar(ncs, 0),
        rho_hat(ncs, 0),
        tau_hat(ncs, 1),
        att(ncs, 1),
        hs(ncs, 0),
        n_chs(ncs, 0)
        { Assign(gt); }
    void Assign(const gen_type& gt,
        flt a_rho0 = 0, flt a_tau0 = 1, flt a_tau = 16, flt a_alph_hs = 0.25);
    // Make choice using Thompson sampling (srn should be vector of ncs
    //  standard normal random numbers)
    void Make_choice(const vf_type& srn);
    // Meta-learning and learning update of agent
    void Update_agent();
    // Update of estimate of recent rewards
    bool Female() const { return female; }
    // public data members
    flt lambda;   // ewm factor for rewards
    flt thr1;     // threshold for turning off attention
    flt thr2;     // threshold for turning on attention
    flt t_thr;    // threshold for number of rounds since attention was off
    flt rho0;     // mean of log-scale prior feature values
    flt tau0;     // precision of log-scale prior feature values
    flt tau;      // precision of log-scale reward variation
    flt alph_hs;  // habit strength (consistency of choice) learning rate
    flt Rew;      // perceived reward
    flt R_ewm;    // estimate of recent rewards
    flt R_tot;    // lifetime reward (so far)
    int ncs;      // number of stimulus dimensions (= number of compound stims)
    int choice;   // choice
    int n_Rew;    // number of rewards so far
    int n_off;    // number of stimulus dimensions with attention off
    int i_num;    // individual number in local group
    int g_num;    // local group number
    bool female;
    bool alive;
    vf_type y_bar;     // mean log-scale rewards
    vf_type rho_hat;   // log-scale mean parameters
    vf_type tau_hat;   // log-scale precision parameters
    vf_type att;       // attention (1/0) to stimulus dimensions
    vf_type hs;        // consistency of choice (habit strength) estimate
    vi_type n_chs;     // number of times CS has been chosen
};

template<typename GenType>
void Phenotype<GenType>::Assign(const gen_type& gt,
    flt a_rho0, flt a_tau0, flt a_tau, flt a_alph_hs)
{
    val_type val = gt.Value();
    // assume val is a vector with components corresponding to the traits 
    // lambda, thr1, thr2, t_thr
    lambda = val[0];
    thr1 = val[1];
    thr2 = val[2];
    t_thr = val[3];
    rho0 = a_rho0;
    tau0 = a_tau0;
    tau = a_tau;
    alph_hs = a_alph_hs;
    Rew = 0;
    R_ewm = 0;
    R_tot = 0;
    choice = 0;
    n_Rew = 0;
    n_off = 0;
    i_num = 0;
    g_num = 0;
    female = false;
    alive = false;
    for (int k = 0; k < ncs; ++k) {
        y_bar[k] = 0;
        rho_hat[k] = rho0;
        tau_hat[k] = tau0;
        att[k] = 1;
        hs[k] = 0;
        n_chs[k] = 0;
    }
}

// Choose option using the TS procedure
template<typename GenType>
void Phenotype<GenType>::Make_choice(const vf_type& srn) {
    // vector of random draws from log-scale posteriors
    flt lo_val = -4/std::sqrt(tau0);
    vf_type y_hat(ncs, lo_val);
    for (int k = 0; k < ncs; ++k) {
        if (att[k] > 0) {
            y_hat[k] = rho_hat[k] + srn[k]/std::sqrt(tau_hat[k]);
        }
    }
    if (n_off < ncs) {
        choice = std::distance(y_hat.begin(), 
            std::max_element(y_hat.begin(), y_hat.end()));
    }
    // if all attention is off, then keep the previous choice
}

// Meta-learning update of agent state. The agent first updates the
// consistency of choice (for the chosen option), then updates the ewm
// estimate of rewards, then updates posterior parameters (TS procedure), then
// decides whether to turn off attention. The agent then checks if recent
// rewards have become so small that all attention should be turned on.
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

    // update posterior parameters
    if (att[choice] > 0) {
        // update mean log-scale reward
        y_bar[choice] += (std::log(Rew) - y_bar[choice])/(n_chs[choice] + 1);
        // update number of times this feature has been chosen
        ++n_chs[choice];
        // update log-scale mean parameter
        rho_hat[choice] = (rho_hat[choice]*tau_hat[choice] + 
            tau*y_bar[choice])/(tau_hat[choice] + tau);
        // update log-scale precision parameter
        tau_hat[choice] += tau;
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
                            ++n_off;
                        }
                    }
                }
            }
        }
    }
    // check if all attention should be turned on (with reset of state)
    if (n_off > 0 && n_Rew > t_thr && R_ewm < thr2) {
        // turn on attention and reset for all dimensions
        for (int k = 0; k < ncs; ++k) {
            // reset mean log-scale rewards
            y_bar[k] = 0;
            // reset log-scale mean and precision parameters
            rho_hat[k] = rho0;
            tau_hat[k] = tau0;
            // set attention to on
            att[k] = 1;
            // reset consistency of choice to zero
            hs[k] = 0;
            // reset count of feature observations
            n_chs[k] = 0;
        }
        // attention is on for all dimensions
        n_off = 0;
        // reset count of rewards
        n_Rew = 0;
    }
    // update total rewards (expected reproductive success)
    R_tot += Rew;
}

#endif // PHENOTYPE_HPP
