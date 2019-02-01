//
//  lpm_likelihood.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <math.h>

#include <spf/numerical_utils.hpp>

#include "lpm_likelihood.hpp"

//double *construct_mutation_matrix(vector<size_t> pathway_membership, size_t n_patients, size_t n_genes, size_t n_pathways)
//{
//    // construct mutation matrix
//    // calculate sizes for each pathways as well
//    double *mut_matrix = new double[n_patients*n_pathways];
//    for (size_t n = 0; n < n_genes; n++) {
//        size_t membership_idx = pathway_membership.at(n);
//        for (size_t k = 0; k < n_pathways; k++) {
//            mut_matrix[n*n_pathways + k] = 0.0;
//            if (membership_idx == k) {
//                mut_matrix[n*n_pathways + k] = 1.0;
//            }
//        }
//    }
//    return mut_matrix;
//}
//
//void fast_matrix_product(gsl_matrix_view &obs, size_t m, vector<size_t> &pathway_membership, size_t num_pathways, vector<size_t> &r)
//{
//    size_t N = pathway_membership.size();
//    for (size_t n = 0; n < N; n++) {
//        double val = gsl_matrix_get(&obs.matrix, m, n);
//        if (val == 1.0) {
//            size_t pathway_idx = pathway_membership[n];
//            if (pathway_idx >= r.size()) {
//                cerr << "Error in fast_matrix_product: invalid size of ret vector is allocated." << endl;
//                exit(-1);
//            }
//            r[pathway_idx] += 1;
//        }
//    }
//}

double compute_log_lik_active(double n_mutations, double pathway_size, double bgp, double fbp)
{
    if (pathway_size <= 0) {
        // this is not in the support set
        return DOUBLE_NEG_INF;
    }
    double log_lik = 0.0;
    if (n_mutations == 0) {
        // one entry got flipped back from 1 -> 0 and the rest are not impacted by back ground prob
        log_lik += log(fbp);
        log_lik += (pathway_size - 1) * log(1-bgp);
        log_lik += log(pathway_size);
    } else if (n_mutations == 1) {
        // no error case: one entry was unaffected by fbp (1->1) and the rest are not affected by background prob
        log_lik += log(1 - fbp);
        log_lik += (pathway_size - 1) * log(1 - bgp);
        if (pathway_size > n_mutations) {
            // we need to account for flip back (1->0) of a gene
            // followed by background mutation of the one that gave rise to n_mutations = 1
            // remaining genes are unaffected
            // e.g., obs may be 1,0,0, which may have been obtained by
            // fb and bg mutation from 0,1,0 or 0,0,1
            double log_lik2 = log(fbp) + log(bgp);
            log_lik2 += (pathway_size - 2) * log(1 - bgp);
            log_lik2 += log(pathway_size - 1);
            log_lik = log_add(log_lik, log_lik2);
        }
    } else {
        // one entry was unaffected by fbp (1->1)
        // (n_mutations - 1) are background mutations
        // the rest are not affected by background prob
        // e.g., 1,1,1,0 can be obtained from 1,0,0,0 with bgm on second and third entry
        // but we do not know which one is the unaffected one, so the number of possibilities are
        // (1,0,0,0), (0,1,0,0) or (0,0,1,0) with 2 bg mutations, no flip back, and no bg on the last entry.
        log_lik += log(n_mutations); // multiply by the number of possibilities
        log_lik += log(1 - fbp);
        log_lik += (n_mutations - 1) * log(bgp);
        log_lik += (pathway_size - n_mutations) * log(1 - bgp);
        if (pathway_size > n_mutations) {
            // account for flip back (1->0) followed by n_mutations bgm
            // e.g., 0,0,0,1->0,0,0,0->1,1,1,0
            // note: we don't allow fb followed by bgm
            double log_lik2 = log(pathway_size - n_mutations);
            log_lik2 += log(fbp);
            log_lik2 += n_mutations * log(bgp);
            log_lik2 += (pathway_size - n_mutations - 1) * log(1 - bgp);
            log_lik = log_add(log_lik, log_lik2);
        }
    }
    log_lik -= log(pathway_size);
    return log_lik;
}

double compute_log_lik_not_active(double n_mutations, double pathway_size, double bgp)
{
    // all mutations are background mutations
    double log_lik = n_mutations * log(bgp);
    log_lik += (pathway_size - n_mutations) * log(1 - bgp);
    return log_lik;
}

double compute_likelihood_for_sample(vector<size_t> &r, LinearProgressionState &state, size_t stage, double bgp, double fbp)
{
    double log_lik = 0.0;
    size_t K = state.get_num_pathways();
    for (int k = 0; k < K; k++) {
        double val = r[k];
        if (k <= stage) {
            log_lik += compute_log_lik_active(val, state.get_pathway_size(k), bgp, fbp);
        } else {
            log_lik += compute_log_lik_not_active(val, state.get_pathway_size(k), bgp);
        }
    }
    return log_lik;
}

void reset_r_vector(vector<size_t> &vec)
{
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = 0;
    }
}

bool contains_empty_pathway(vector<size_t> &pathway_sizes)
{
    for (size_t i = 0; i < pathway_sizes.size(); i++) {
        if (pathway_sizes[i] == 0) {
            return true;
        }
    }
    return false;
}

double compute_pathway_likelihood(gsl_matrix_view &obs, vector<size_t> &sum_vec, LinearProgressionState &state, LinearProgressionParameters &params)
{
    double log_lik = 0.0;

    unsigned int M = obs.matrix.size1; // number of observations
    //unsigned int N = obs.matrix.size2; // number of genes
    unsigned int K = state.get_num_pathways(); // number of driver pathways

    if (state.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }

    vector<size_t> *stages = params.has_patient_stages() ? &params.get_patient_progression_stages() : 0;
    vector<size_t> ret(K);
    if (stages != 0) {
        for (size_t m = 0; m < M; m++) {
            //fast_matrix_product(obs, m, pathway_membership, K, ret);
            state.compute_counts_for_sample(obs, sum_vec, m, ret);
            double log_lik_m = compute_likelihood_for_sample(ret, state, stages->at(m), params.get_bgp(), params.get_fbp());
            log_lik += log_lik_m;
        }
    } else {
        // marginalize out stages for each patient
        for (size_t m = 0; m < M; m++) {
            double log_lik_m = DOUBLE_NEG_INF;
            //fast_matrix_product(obs, m, pathway_membership, K, ret);
            state.compute_counts_for_sample(obs, sum_vec, m, ret);
            // marginalize over patient stages taking on {1, ..., driver pathways}
            for (size_t stage = 0; stage < state.get_num_driver_pathways(); stage++) {
                double log_lik_stage = compute_likelihood_for_sample(ret, state, stage, params.get_bgp(), params.get_fbp());
                log_lik_stage -= log(K);
                log_lik_m = log_add(log_lik_m, log_lik_stage);
            }
            log_lik += log_lik_m;
        }
    }
    return log_lik;
}

double compute_pathway_likelihood(vector<vector<size_t>> &R, LinearProgressionState &state, LinearProgressionParameters &params)
{
    if (R.size() == 0) {
        cerr << "Error: R matrix is empty" << endl;
        exit(-1);
    }
    
    double log_lik = 0.0;
    unsigned int M = R.size(); // number of observations
    unsigned int K = state.get_num_pathways(); // number of driver pathways
    
    if (state.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }
    
    vector<size_t> *stages = params.has_patient_stages() ? &params.get_patient_progression_stages() : 0;
    if (stages != 0) {
        for (size_t m = 0; m < M; m++) {
            double log_lik_m = compute_likelihood_for_sample(R[m], state, stages->at(m), params.get_bgp(), params.get_fbp());
            log_lik += log_lik_m;
        }
    } else {
        // marginalize out stages for each patient
        for (size_t m = 0; m < M; m++) {
            double log_lik_m = DOUBLE_NEG_INF;
            // marginalize over patient stages taking on {1, ..., driver pathways}
            for (size_t stage = 0; stage < state.get_num_driver_pathways(); stage++) {
                double log_lik_stage = compute_likelihood_for_sample(R[m], state, stage, params.get_bgp(), params.get_fbp());
                log_lik_stage -= log(K);
                log_lik_m = log_add(log_lik_m, log_lik_stage);
            }
            log_lik += log_lik_m;
        }
    }
    return log_lik;
}



