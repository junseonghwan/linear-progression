//
//  lpm_pg_proposal.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-27.
//

#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_pg_proposal.hpp"

LPMParamProposal::LPMParamProposal(size_t n_mh_iters, double fbp_max, double bgp_max, double mh_proposal_sd) :
n_mh_iter(n_mh_iters), fbp_max(fbp_max), bgp_max(bgp_max), mh_proposal_sd_fbp(mh_proposal_sd), mh_proposal_sd_bgp(mh_proposal_sd)
{
    
}

shared_ptr<LinearProgressionParameters> LPMParamProposal::sample_from_prior(gsl_rng *random)
{
    double bgp = gsl_ran_flat(random, 0.0, bgp_max);
    double fbp = bgp;
    if (fbp_max > 0.0) {
        fbp = gsl_ran_flat(random, 0.0, fbp_max);
    }
    cout << "Initial params: fbp=" << fbp << ", bgp=" << bgp << endl;
    shared_ptr<LinearProgressionParameters> ret(new LinearProgressionParameters(fbp, bgp));
    return ret;
}

void LPMParamProposal::sample_separately(gsl_rng *random, const LinearProgressionState &state, LinearProgressionParameters &new_param)
{
    if (bgps.size() < 100) {
        mh_proposal_sd_bgp = new_param.get_bgp()/2;
        mh_proposal_sd_fbp = new_param.get_fbp()/2;
    }

    double old_log_lik = compute_pathway_likelihood(state, new_param);
    double old_bgp = new_param.get_bgp();
    double old_fbp = new_param.get_fbp();
    cout << "(" << old_fbp << ", " << old_bgp << ") -> ";

    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_bgp = gsl_ran_gaussian(random, mh_proposal_sd_bgp) + old_bgp;
        if (new_bgp > 0 && new_bgp <= bgp_max) {
            new_param.set_bgp(new_bgp);
            double new_log_lik = compute_pathway_likelihood(state, new_param);
            double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
            if (log_u < (new_log_lik - old_log_lik)) {
                // accept
                old_log_lik = new_log_lik;
                old_bgp = new_bgp;
            } else {
                // revert
                new_param.set_bgp(old_bgp);
            }
        }
        bgps.push_back(old_bgp);
    }
    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_fbp = gsl_ran_gaussian(random, mh_proposal_sd_fbp) + old_fbp;
        if (new_fbp > 0 && new_fbp <= fbp_max) {
            new_param.set_fbp(new_fbp);
            double new_log_lik = compute_pathway_likelihood(state, new_param);
            double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
            if (log_u < (new_log_lik - old_log_lik)) {
                // accept
                old_log_lik = new_log_lik;
                old_fbp = new_fbp;
            } else {
                // revert
                new_param.set_fbp(old_fbp);
            }
        }
        fbps.push_back(old_fbp);
    }
    cout << "(" << new_param.get_fbp() << ", " << new_param.get_bgp() << ")" << endl;
}

void LPMParamProposal::sample_together(gsl_rng *random, const LinearProgressionState &state, LinearProgressionParameters &new_param)
{
    if (bgps.size() < 100) {
        mh_proposal_sd_bgp = new_param.get_bgp()/2;
    }

    double old_log_lik = compute_pathway_likelihood(state, new_param);
    double old_error_prob = new_param.get_bgp();
    // check:
    if (new_param.get_bgp() != new_param.get_fbp()) {
        cerr << "Error: 1 parameter model but FBP != BGP: " << new_param.get_fbp() << " != " << new_param.get_bgp() << endl;
        exit(-1);
    }
    cout << "Curr state: " << state.to_string() << endl;

    for (size_t i = 0; i < n_mh_iter; i++) {
        cout << "(" << new_param.get_fbp() << ", " << old_log_lik << ")" << endl;
        double new_error_prob = gsl_ran_gaussian(random, mh_proposal_sd_bgp) + old_error_prob;
        if (new_error_prob > 0 && new_error_prob <= bgp_max) {
            new_param.set_bgp(new_error_prob);
            new_param.set_fbp(new_error_prob);
            double new_log_lik = compute_pathway_likelihood(state, new_param);
            double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
            if (log_u < (new_log_lik - old_log_lik)) {
                // accept
                old_log_lik = new_log_lik;
                old_error_prob = new_error_prob;
            } else {
                // revert
                new_param.set_bgp(old_error_prob);
                new_param.set_fbp(old_error_prob);
            }
        }
        bgps.push_back(old_error_prob);
        fbps.push_back(old_error_prob);
    }

    cout << "(" << new_param.get_bgp() << ", " << old_log_lik << ")" << endl;
}

shared_ptr<LinearProgressionParameters> LPMParamProposal::propose(gsl_rng *random, const LinearProgressionParameters &curr, shared_ptr<ParticleGenealogy<LinearProgressionState>> genealogy)
{
    LinearProgressionParameters *new_param = new LinearProgressionParameters(curr.get_fbp(), curr.get_bgp());

    const LinearProgressionState &state = genealogy->get_state_at(genealogy->size()-1);

    // 1. fbp_max = 0.0 => fbp = bgp
    // 2. fbp_max > 0.0 => fbp, bgp should be sampled separately
    if (fbp_max == 0.0) {
        sample_together(random, state, *new_param);
    } else {
        sample_separately(random, state, *new_param);
    }

    return shared_ptr<LinearProgressionParameters>(new_param);
}

double LPMParamProposal::log_prior(const LinearProgressionParameters &curr)
{
    // proportional to 0 (Uniform prior with support (0, fbp_max], (0, bgp_max]).
    return 0.0;
}
