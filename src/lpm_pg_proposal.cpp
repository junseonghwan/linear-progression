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

LPMParamProposal::LPMParamProposal(gsl_matrix_view &obs_view, vector<size_t> &row_sum, size_t n_patients, size_t n_pathways, size_t n_genes, size_t n_mh_iters, double a, double b) :
n_patients(n_patients), n_pathways(n_pathways), n_genes(n_genes), n_mh_iter(n_mh_iters), obs_view(obs_view), row_sum(row_sum), a(a), b(b)
{
    stages = new vector<size_t>(n_patients, 0);
}

LinearProgressionParameters *LPMParamProposal::sample_from_prior(gsl_rng *random)
{
    double fbp = gsl_ran_flat(random, a, b);
    double bgp = gsl_ran_flat(random, a, b);
    //LinearProgressionParameters *new_param = new LinearProgressionParameters(fbp, bgp, *stages);
    LinearProgressionParameters *new_param = new LinearProgressionParameters(fbp, bgp);
    return new_param;
}

LinearProgressionParameters *LPMParamProposal::propose(gsl_rng *random, LinearProgressionParameters *curr, ParticleGenealogy<LinearProgressionState> *genealogy)
{
    LinearProgressionState &state = genealogy->at(genealogy->size()-1)->first;

    // perform MH on each of the parameters
    double old_fbp = curr->get_fbp();
    double old_bgp = curr->get_bgp();

    // compute the counts using obs, row_sum, and curr
    size_t K = state.get_num_pathways();
    vector<vector<size_t>> R(n_patients);
    for (size_t m = 0; m < n_patients; m++) {
        vector<size_t> r(K);
        state.compute_counts_for_sample(obs_view, row_sum, m, r);
        R[m] = r;
    }

    double old_log_lik = compute_pathway_likelihood(R, state, *curr);
    cout << "(" << old_fbp << ", " << old_bgp << ") -> ";
    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_fbp = gsl_ran_gaussian(random, mh_proposal_sd) + old_fbp;
        if (new_fbp < 0.03) {
            continue;
        }
        
        curr->set_fbp(new_fbp);
        double new_log_lik = compute_pathway_likelihood(R, state, *curr);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - old_log_lik)) {
            // accept
            old_log_lik = new_log_lik;
            old_fbp = new_fbp;
        } else {
            // revert
            curr->set_fbp(old_fbp);
        }
    }
    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_bgp = gsl_ran_gaussian(random, mh_proposal_sd) + old_bgp;
        if (new_bgp < 0.03) {
            continue;
        }

        curr->set_bgp(new_bgp);
        double new_log_lik = compute_pathway_likelihood(R, state, *curr);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - old_log_lik)) {
            // accept
            old_log_lik = new_log_lik;
            old_bgp = new_bgp;
        } else {
            // revert
            curr->set_bgp(old_bgp);
        }
    }
    cout << "(" << curr->get_fbp() << ", " << curr->get_bgp() << ")" << endl;
    return curr;
}

double LPMParamProposal::log_prior(LinearProgressionParameters *curr)
{
    // proportional to 0
    return 0.0;
}
