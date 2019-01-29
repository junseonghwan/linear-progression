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

// initialize and maintain the vector over the patient stages
//vector<size_t> *stages = new vector<size_t>(n_patients);
//for (unsigned int i = 0; i < n_patients; i++) {
//    unsigned int k = gsl_rng_uniform_int(random, n_pathways);
//    (*stages)[k] = i;
//}

LPMParamProposal::LPMParamProposal(gsl_matrix_view &obs_view, size_t n_patients, size_t n_pathways, size_t n_genes, size_t n_mh_iters, double a, double b) :
n_patients(n_patients), n_pathways(n_pathways), n_genes(n_genes), n_mh_iter(n_mh_iters), obs_view(obs_view), a(a), b(b)
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

//    // compute matrix product
//    vector<size_t> pathway_membership = state.get_pathway();
//    double mut_matrix[n_genes*n_pathways];
//    vector<size_t> pathway_sizes(n_pathways, 0);
//    for (size_t n = 0; n < n_patients; n++) {
//        size_t membership_idx = pathway_membership.at(n);
//        pathway_sizes[membership_idx] += 1;
//        for (size_t k = 0; k < n_pathways; k++) {
//            mut_matrix[n*n_pathways + k] = 0.0;
//            if (membership_idx == k) {
//                mut_matrix[n*n_pathways + k] = 1.0;
//            }
//        }
//    }
//    gsl_matrix_view mut_matrix_view = gsl_matrix_view_array(mut_matrix, n_genes, n_pathways);
//    // compute matrix multiplication
//    double ret_matrix[n_patients*n_pathways];
//    gsl_matrix_view ret_matrix_view = gsl_matrix_view_array(ret_matrix, n_patients, n_pathways);
//    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
//                    1.0, &obs_view.matrix, &mut_matrix_view.matrix,
//                    0.0, &ret_matrix_view.matrix);
//
//    for (size_t n = 0; n < n_genes; n++) {
//        size_t membership_idx = pathway_membership.at(n);
//        pathway_sizes[membership_idx] += 1;
//    }
//    // sample patient stages
//    vector<double> log_liks(n_pathways);
//    vector<double> probs(n_pathways);
//    for (size_t m = 0; m < n_patients; m++) {
//        for (size_t k = 0; k < n_pathways; k++) {
//            log_liks[k] = compute_likelihood(ret_matrix_view, m, k, pathway_sizes, curr->get_bgp(), curr->get_fbp());
//        }
//        normalize(log_liks, probs);
//        size_t kstar = multinomial(random, probs);
//        (*stages)[m] = kstar;
//    }

    // perform MH on each of the parameters
    double old_fbp = curr->get_fbp();
    double old_bgp = curr->get_bgp();
    double old_log_lik = compute_pathway_likelihood(obs_view, state, *curr);
    cout << "(" << old_fbp << ", " << old_bgp << ") -> ";
    for (size_t i = 0; i < n_mh_iter; i++) {
        double new_fbp = gsl_ran_flat(random, a, b);
        double new_bgp = gsl_ran_flat(random, a, b);
        curr->set_fbp(new_fbp);
        curr->set_bgp(new_bgp);
        double new_log_lik = compute_pathway_likelihood(obs_view, state, *curr);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - old_log_lik)) {
            // accept
            old_log_lik = new_log_lik;
            old_fbp = new_fbp;
            old_bgp = new_bgp;
        } else {
            // revert -- nothing to do
        }
    }
    cout << "(" << old_fbp << ", " << old_bgp << ")" << endl;
    curr->set_fbp(old_fbp);
    curr->set_bgp(old_bgp);
    return curr;
}

double LPMParamProposal::log_prior(LinearProgressionParameters *curr)
{
    // proportional to 0
    return 0.0;
}
