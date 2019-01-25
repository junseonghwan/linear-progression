//
//  lpm_likelihood.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <math.h>

#include <spf/numerical_utils.hpp>

#include "lpm_likelihood.hpp"

double compute_log_lik(double n_mutations, double pathway_size, double bgp, double fbp)
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

double compute_log_lik_bg(double n_mutations, double pathway_size, double bgp, double fbp)
{
    // all mutations are background mutations
    double log_lik = n_mutations * log(bgp);
    log_lik += (pathway_size - n_mutations) * log(1 - bgp);
    return log_lik;
}

double compute_pathway_likelihood(gsl_matrix_view &Y, gsl_matrix_view &X, vector<size_t> &stages, double fbp, double bgp)
{
    double log_lik = 0.0;

    unsigned int M = Y.matrix.size1; // number of observations
    unsigned int N = Y.matrix.size2; // number of genes
    unsigned int K = X.matrix.size2; // number of pathways
    
    // count pathway sizes
    vector<int> pathway_sizes(K, 0);
    for (size_t i = 0; i < N; i++) {
        for (size_t k = 0; k < K; k++) {
            if (gsl_matrix_get(&X.matrix, i, k) == 1) {
                pathway_sizes[k] += 1;
                break;
            }
        }
    }

    double r[M*K];
    gsl_matrix_view R = gsl_matrix_view_array(r, M, K);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, &Y.matrix, &X.matrix,
                    0.0, &R.matrix);

    for (int m = 0; m < M; m++) {
        double log_lik_m = 0.0;
        for (int k = 0; k < K; k++) {
            double val = gsl_matrix_get(&R.matrix, m, k);
            if (k <= stages[m]) {
                log_lik_m += compute_log_lik(val, pathway_sizes[k], bgp, fbp);
            } else {
                log_lik_m += compute_log_lik_bg(val, pathway_sizes[k], bgp, fbp);
            }
        }
        log_lik += log_lik_m;
    }
    
    return log_lik;
}

double compute_likelihood(gsl_matrix_view &ret, size_t m, size_t stage, vector<size_t> pathway_sizes, double bgp, double fbp)
{
    double log_lik_m = 0.0;
    size_t K = pathway_sizes.size();
    for (int k = 0; k < K; k++) {
        double val = gsl_matrix_get(&ret.matrix, m, k);
        if (k <= stage) {
            log_lik_m += compute_log_lik(val, pathway_sizes[k], bgp, fbp);
        } else {
            log_lik_m += compute_log_lik_bg(val, pathway_sizes[k], bgp, fbp);
        }
    }
    return log_lik_m;
}

double compute_pathway_likelihood(gsl_matrix_view &obs, LinearProgressionState &state, LinearProgressionParameters &params)
{
    double log_lik = 0.0;

    unsigned int M = obs.matrix.size1; // number of observations
    unsigned int N = obs.matrix.size2; // number of genes
    unsigned int K = state.get_num_pathways(); // number of pathways

    // construct mutation matrix
    // calculate sizes for each pathways as well
    vector<size_t> pathway_membership = state.get_pathway();
    double mut_matrix[N*K];
    vector<size_t> pathway_sizes(K, 0);

    for (size_t n = 0; n < N; n++) {
        size_t membership_idx = pathway_membership.at(n);
        if (membership_idx >= N) {
            cerr << "Error" << endl;
            exit(-1);
        }
        pathway_sizes[membership_idx] += 1;
        for (size_t k = 0; k < K; k++) {
            mut_matrix[n*K + k] = 0.0;
            if (membership_idx == k) {
                mut_matrix[n*K + k] = 1.0;
            }
        }
    }
    
    // if pathway size == 0 for any, return -inf
    for (size_t k = 0; k < K; k++) {
        if (pathway_sizes[k] == 0)
            return DOUBLE_NEG_INF;
    }
    
    gsl_matrix_view mut_matrix_view = gsl_matrix_view_array(mut_matrix, N, K);

    // compute matrix multiplication
    double ret_matrix[M*K];
    gsl_matrix_view ret_matrix_view = gsl_matrix_view_array(ret_matrix, M, K);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, &obs.matrix, &mut_matrix_view.matrix,
                    0.0, &ret_matrix_view.matrix);

    // using the patient stages available in the params, complete the likelihood calculation
    vector<size_t> *stages = 0;
    if (params.has_patient_stages()) {
         stages = &params.get_patient_progression_stages();
    }
    if (stages != 0) {
        for (int m = 0; m < M; m++) {
            double log_lik_m = 0.0;
            for (int k = 0; k < K; k++) {
                double val = gsl_matrix_get(&ret_matrix_view.matrix, m, k);
                if (k <= stages->at(m)) {
                    log_lik_m += compute_log_lik(val, pathway_sizes[k], params.get_bgp(), params.get_fbp());
                } else {
                    log_lik_m += compute_log_lik_bg(val, pathway_sizes[k], params.get_bgp(), params.get_fbp());
                }
            }
            log_lik += log_lik_m;
        }
    } else {
        // marginalize out stages for each patient
        for (size_t m = 0; m < M; m++) {
            double log_lik_m = DOUBLE_NEG_INF;
            for (size_t stage = 0; stage < K; stage++) {
                double log_lik_stage = compute_likelihood(ret_matrix_view, m, stage, pathway_sizes, params.get_bgp(), params.get_fbp());
                log_lik_m = log_add(log_lik_m, log_lik_stage);
            }
            log_lik += log_lik_m;
        }
    }

    return log_lik;
}
