//
//  lpm_pg_proposal.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-27.
//

#include <iostream>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_pg_proposal.hpp"

LPMParamProposal::LPMParamProposal(gsl_matrix &obs_matrix, vector<size_t> &row_sum, size_t n_patients, size_t n_mh_iters, double fbp_max, double bgp_max) :
n_patients(n_patients), n_mh_iter(n_mh_iters), obs_matrix(obs_matrix), row_sum(row_sum), fbp_max(fbp_max), bgp_max(bgp_max)
{
    stages = new vector<size_t>(n_patients, 0);
}

LinearProgressionParameters *LPMParamProposal::sample_from_prior(gsl_rng *random)
{
    double fbp = gsl_ran_flat(random, 0.0, fbp_max);
    double bgp = gsl_ran_flat(random, 0.0, bgp_max);
    //LinearProgressionParameters *new_param = new LinearProgressionParameters(fbp, bgp, *stages);
    LinearProgressionParameters *new_param = new LinearProgressionParameters(fbp, bgp);
    return new_param;
}

LinearProgressionParameters *LPMParamProposal::propose(gsl_rng *random, LinearProgressionParameters *curr, ParticleGenealogy<LinearProgressionState> *genealogy)
{
    LinearProgressionParameters *new_param = new LinearProgressionParameters(curr->get_fbp(), curr->get_bgp());

    const LinearProgressionState &state = genealogy->get_state_at(genealogy->size()-1);

    // perform MH on each of the parameters
    double old_fbp = new_param->get_fbp();
    double old_bgp = new_param->get_bgp();

    // compute the counts using obs, row_sum, and curr
    size_t K = state.get_num_pathways();
    vector<vector<size_t>> R(n_patients);
    for (size_t m = 0; m < n_patients; m++) {
        vector<size_t> r(K);
        state.compute_counts_for_sample(obs_matrix, row_sum, m, r);
        R[m] = r;
    }

    double old_log_lik = compute_pathway_likelihood(R, state, *curr);
    cout << "(" << old_fbp << ", " << old_bgp << ") -> ";
    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_fbp = gsl_ran_gaussian(random, mh_proposal_sd) + old_fbp;
        if (new_fbp < 0 || new_fbp > fbp_max) {
            continue;
        }
        
        new_param->set_fbp(new_fbp);
        double new_log_lik = compute_pathway_likelihood(R, state, *curr);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - old_log_lik)) {
            // accept
            old_log_lik = new_log_lik;
            old_fbp = new_fbp;
        } else {
            // revert
            new_param->set_fbp(old_fbp);
        }
    }
    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_bgp = gsl_ran_gaussian(random, mh_proposal_sd) + old_bgp;
        if (new_bgp < 0 || new_bgp > bgp_max) {
            continue;
        }

        new_param->set_bgp(new_bgp);
        double new_log_lik = compute_pathway_likelihood(R, state, *curr);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - old_log_lik)) {
            // accept
            old_log_lik = new_log_lik;
            old_bgp = new_bgp;
        } else {
            // revert
            new_param->set_bgp(old_bgp);
        }
    }
    cout << "(" << new_param->get_fbp() << ", " << new_param->get_bgp() << ")" << endl;
    return new_param;
}

double LPMParamProposal::log_prior(LinearProgressionParameters *curr)
{
    // proportional to 0
    return 0.0;
}
