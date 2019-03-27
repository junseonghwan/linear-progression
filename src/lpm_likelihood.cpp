//
//  lpm_likelihood.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <math.h>

#include <spf/numerical_utils.hpp>

#include "lpm_likelihood.hpp"

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

double compute_log_lik_inactive(double n_mutations, double pathway_size, double bgp)
{
    // all mutations are background mutations
    double log_lik = n_mutations * log(bgp);
    log_lik += (pathway_size - n_mutations) * log(1 - bgp);
    return log_lik;
}

double compute_likelihood_for_sample(const vector<size_t> &r, const LinearProgressionState &state, size_t stage, double bgp, double fbp)
{
    double log_lik = 0.0;
    size_t K = state.get_num_pathways();
    for (int k = 0; k < K; k++) {
        double val = r[k];
        if (k <= stage) {
            log_lik += compute_log_lik_active(val, state.get_pathway_size(k), bgp, fbp);
        } else {
            log_lik += compute_log_lik_inactive(val, state.get_pathway_size(k), bgp);
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

//double compute_pathway_likelihood(gsl_matrix &obs_matrix,
//                                  vector<size_t> &sum_vec,
//                                  const LinearProgressionState &state,
//                                  const LinearProgressionParameters &params,
//                                  unordered_map<size_t, unordered_map<size_t, double> > _dict_active_probs,
//                                  unordered_map<size_t, unordered_map<size_t, double> > _dict_inactive_probs)
//{
//    double log_lik = 0.0;
//
//    unsigned int n_obs = obs_matrix.size1; // number of observations
//    //unsigned int N = obs.matrix.size2; // number of genes
//    size_t n_pathways = state.get_num_pathways();
//
//    if (state.contains_empty_driver_pathway()) {
//        return DOUBLE_NEG_INF;
//    }
//    double log_n_pathways = log(n_pathways);
//
//    for (size_t m = 0; m < n_obs; m++) {
//        const vector<size_t> &ret = state.get_cache_at(m);
//        double log_lik_m = 0.0;
//
//        // marginalize out the patient progression stages
//        // first, compute the log_lik for stage = 0
//        double log_lik_stage = compute_likelihood_for_sample(ret, state, 0, params.get_bgp(), params.get_fbp());
//        log_lik_m = log_lik_stage - log_n_pathways; // -log_n_pathways accounts for prior on stage assignment: (1/n_pathways)
//
//        // second, compute the log_lik for stage 1 up to number of driver pathways
//        for (size_t stage = 1; stage < state.get_num_driver_pathways(); stage++) {
//            size_t r = ret[stage];
//            size_t pathway_size = state.get_pathway_size(stage);
//            if (!_dict_inactive_probs[pathway_size].count(r)) {
//                _dict_inactive_probs[pathway_size][r] = compute_log_lik_inactive(r, pathway_size, params.get_bgp());
//            }
//            if (!_dict_active_probs[pathway_size].count(r)) {
//                _dict_active_probs[pathway_size][r] = compute_log_lik_active(r, pathway_size, params.get_bgp(), params.get_fbp());
//            }
//            log_lik_stage -= _dict_inactive_probs.at(pathway_size).at(r);
//            log_lik_stage += _dict_active_probs.at(pathway_size).at(r);
//            log_lik_m = log_add(log_lik_m, log_lik_stage - log_n_pathways);
//        }
//        log_lik += log_lik_m;
//    }
//
//    return log_lik;
//}

double compute_pathway_likelihood(gsl_matrix &obs_matrix,
                                  vector<size_t> &sum_vec,
                                  const LinearProgressionState &state,
                                  const LinearProgressionParameters &params)
{
    double log_lik = 0.0;

    unsigned int n_obs = obs_matrix.size1; // number of observations
    //unsigned int N = obs.matrix.size2; // number of genes
    unsigned int n_pathways = state.get_num_pathways(); // number of driver pathways
    double log_n_pathways = log(n_pathways);
    if (state.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }
    
    //marginalize out stages for each patient
    for (size_t m = 0; m < n_obs; m++) {
        const vector<size_t> &ret = state.get_cache_at(m);
        double log_lik_stage = compute_likelihood_for_sample(ret, state, 0, params.get_bgp(), params.get_fbp());
        double log_lik_m = log_lik_stage - log_n_pathways; // -log_n_pathways accounts for prior on stage assignment: (1/n_pathways)
        // marginalize over patient stages taking on {1, ..., driver pathways}
        for (size_t stage = 1; stage < state.get_num_driver_pathways(); stage++) {
            size_t pathway_size = state.get_pathway_size(stage);
            size_t r = ret[stage];
            log_lik_stage += compute_log_lik_active(r, pathway_size, params.get_bgp(), params.get_fbp());
            log_lik_stage -= compute_log_lik_inactive(r, pathway_size, params.get_bgp());
            log_lik_m = log_add(log_lik_m, log_lik_stage - log_n_pathways);
        }
        log_lik += log_lik_m;
    }
    return log_lik;
}

double compute_pathway_likelihood(vector<vector<size_t> > &R, const LinearProgressionState &state, const LinearProgressionParameters &params)
{
    if (R.size() == 0) {
        cerr << "Error: R matrix is empty" << endl;
        exit(-1);
    }
    
    double log_lik = 0.0;
    unsigned int n_obs = R.size(); // number of observations
    unsigned int n_pathways = state.get_num_pathways(); // number of driver pathways
    double log_n_pathways = log(n_pathways);
    
    if (state.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }
    
    // marginalize out stages for each patient
    for (size_t m = 0; m < n_obs; m++) {
        double log_lik_m = DOUBLE_NEG_INF;
        // marginalize over patient stages taking on {1, ..., driver pathways}
        double log_lik_stage = compute_likelihood_for_sample(R[m], state, 0, params.get_bgp(), params.get_fbp());
        for (size_t stage = 1; stage < state.get_num_driver_pathways(); stage++) {
            size_t r = R[m][stage];
            log_lik_stage -= compute_log_lik_inactive(r, state.get_pathway_size(stage), params.get_bgp());
            log_lik_stage += compute_log_lik_active(r, state.get_pathway_size(stage), params.get_bgp(), params.get_fbp());
            log_lik_m = log_add(log_lik_m, log_lik_stage - log_n_pathways);
        }
        log_lik += log_lik_m;
    }
    return log_lik;
}

void fast_matrix_product(gsl_matrix &obs, size_t m, const vector<size_t> &pathway_membership, vector<size_t> &r)
{
    size_t n_genes = pathway_membership.size();
    for (size_t n = 0; n < n_genes; n++) {
        double val = gsl_matrix_get(&obs, m, n);
        if (val == 1.0) {
            size_t pathway_idx = pathway_membership[n];
            if (pathway_idx >= r.size()) {
                cerr << "Error in fast_matrix_product: invalid size of ret vector is allocated." << endl;
                exit(-1);
            }
            r[pathway_idx] += 1;
        }
    }
}


//double compute_pathway_likelihood(gsl_matrix &obs_matrix, vector<size_t> &sum_vec, const LinearProgressionState &state, LinearProgressionParameters &params)
//{
//    double log_lik = 0.0;
//
//    unsigned int M = obs_matrix.size1; // number of observations
//    //unsigned int N = obs.matrix.size2; // number of genes
//    unsigned int K = state.get_num_pathways(); // number of driver pathways
//
//    if (state.contains_empty_driver_pathway()) {
//        return DOUBLE_NEG_INF;
//    }
//
//    vector<size_t> *stages = params.has_patient_stages() ? &params.get_patient_progression_stages() : 0;
//    vector<size_t> ret(K);
//    if (stages != 0) {
//        for (size_t m = 0; m < M; m++) {
//            //fast_matrix_product(obs, m, pathway_membership, K, ret);
//            state.compute_counts_for_sample(obs_matrix, sum_vec, m, ret);
//            double log_lik_m = compute_likelihood_for_sample(ret, state, stages->at(m), params.get_bgp(), params.get_fbp());
//            log_lik += log_lik_m;
//        }
//    } else {
//        // marginalize out stages for each patient
//        for (size_t m = 0; m < M; m++) {
//            double log_lik_m = DOUBLE_NEG_INF;
//            //fast_matrix_product(obs, m, pathway_membership, K, ret);
//            state.compute_counts_for_sample(obs_matrix, sum_vec, m, ret);
//            // marginalize over patient stages taking on {1, ..., driver pathways}
//            for (size_t stage = 0; stage < state.get_num_driver_pathways(); stage++) {
//                double log_lik_stage = compute_likelihood_for_sample(ret, state, stage, params.get_bgp(), params.get_fbp());
//                log_lik_stage -= log(K);
//                log_lik_m = log_add(log_lik_m, log_lik_stage);
//            }
//            log_lik += log_lik_m;
//        }
//    }
//    return log_lik;
//}

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

//double compute_pathway_likelihood_cache_version1(gsl_matrix &obs_matrix, vector<size_t> &sum_vec, const LinearProgressionState &state, const LinearProgressionParameters &params)
//{
//    double log_lik = 0.0;
//
//    unsigned int n_obs = obs_matrix.size1; // number of observations
//    //unsigned int N = obs.matrix.size2; // number of genes
//    size_t n_pathways = state.get_num_pathways();
//
//    if (state.contains_empty_driver_pathway()) {
//        return DOUBLE_NEG_INF;
//    }
//
//    // construct likelihood table for the state:
//    // for each pathway, store a map : # mutations -> likelihood
//    vector<unordered_map<size_t, double>> _dict_active_probs(n_pathways);
//    vector<unordered_map<size_t, double>> _dict_inactive_probs(n_pathways);
//
//    double log_n_pathways = log(n_pathways);
//
//    for (size_t m = 0; m < n_obs; m++) {
//        const vector<size_t> &ret = state.get_cache_at(m);
//
//        // Fill the dictionary:
//        // Note: move this out from here, and do it as part of SMC model. incur a cost of O(n_genes^2)
//        double log_lik_stage = 0.0;
//        for (size_t k = 0; k < state.get_num_pathways(); k++) {
//            size_t r = ret[k];
//            if (_dict_active_probs[k].count(r) == 0) {
//                _dict_active_probs[k][r] = compute_log_lik_active(r, state.get_pathway_size(k), params.get_bgp(), params.get_fbp());
//            }
//            if (_dict_inactive_probs[k].count(r) == 0) {
//                _dict_inactive_probs[k][r] = compute_log_lik_inactive(r, state.get_pathway_size(k), params.get_bgp());
//            }
//            if (k == 0)
//                log_lik_stage += _dict_active_probs[k][r];
//            else
//                log_lik_stage += _dict_inactive_probs[k][r];
//        }
//
//        double log_lik_m = log_lik_stage - log_n_pathways; // log_n_pathways accounts for prior on stage assignment
//        for (size_t stage = 1; stage < state.get_num_driver_pathways(); stage++) {
//            size_t r = ret[stage];
//            log_lik_stage -= _dict_inactive_probs[stage][r];
//            log_lik_stage += _dict_active_probs[stage][r];
//            log_lik_m = log_add(log_lik_m, log_lik_stage - log_n_pathways);
//        }
//        log_lik += log_lik_m;
//    }
//
//    return log_lik;
//}
//
//double compute_pathway_likelihood_cache_version2(gsl_matrix &obs_matrix, vector<size_t> &sum_vec, const LinearProgressionState &state, const LinearProgressionParameters &params, unordered_map<size_t, unordered_map<size_t, double>> &_dict_active_probs, unordered_map<size_t, unordered_map<size_t, double>> &_dict_inactive_probs)
//{
//    double log_lik = 0.0;
//
//    unsigned int M = obs_matrix.size1; // number of observations
//    //unsigned int N = obs.matrix.size2; // number of genes
//    unsigned int K = state.get_num_pathways(); // number of driver pathways
//
//    if (state.contains_empty_driver_pathway()) {
//        return DOUBLE_NEG_INF;
//    }
//
//    double log_K = log(K);
//
//    for (size_t m = 0; m < M; m++) {
//        const vector<size_t> &ret = state.get_cache_at(m);
//
//        // fill the dictionary
//        double log_lik_stage = 0.0;
//        for (size_t k = 0; k < state.get_num_pathways(); k++) {
//            size_t r = ret[k];
//            size_t pw_size = state.get_pathway_size(k);
//            if (_dict_active_probs[pw_size].count(r) == 0) {
//                _dict_active_probs[pw_size][r] = compute_log_lik_active(r, pw_size, params.get_bgp(), params.get_fbp());
//            }
//            if (_dict_inactive_probs[pw_size].count(r) == 0) {
//                _dict_inactive_probs[pw_size][r] = compute_log_lik_inactive(r, pw_size, params.get_bgp());
//            }
//            if (k == 0)
//                log_lik_stage += _dict_active_probs[pw_size][r];
//            else
//                log_lik_stage += _dict_inactive_probs[pw_size][r];
//        }
//
//        double log_lik_m = log_lik_stage - log_K; // log_K accounts for prior on stage assignment
//        for (size_t stage = 1; stage < state.get_num_driver_pathways(); stage++) {
//            size_t r = ret[stage];
//            size_t pw_size = state.get_pathway_size(stage);
//            log_lik_stage -= _dict_inactive_probs[pw_size][r];
//            log_lik_stage += _dict_active_probs[pw_size][r];
//            log_lik_m = log_add(log_lik_m, log_lik_stage - log_K);
//        }
//        log_lik += log_lik_m;
//    }
//
//    return log_lik;
//}



