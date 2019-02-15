//
//  lpm_model_cache.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-15.
//

#include "lpm_model_cache.hpp"

//
//  lpm_model.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <math.h>
#include <numeric>
#include <iostream>
#include <list>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "lpm_model_cache.hpp"

LinearProgressionModelCache::LinearProgressionModelCache(size_t num_genes,
                                               size_t num_driver_pathways,
                                               size_t num_smc_iter,
                                               size_t num_mcmc_iter,
                                               gsl_matrix &obs,
                                               vector<size_t> &row_sum,
                                               double pathway_swap_prob,
                                               bool allocate_passenger_pathway,
                                               MoveType move_type) :
num_genes(num_genes), num_driver_pathways(num_driver_pathways), num_smc_iter(num_smc_iter), num_mcmc_iter(num_mcmc_iter),
obs(obs), row_sum(row_sum), pathway_swap_prob(pathway_swap_prob), allocate_passenger_pathway(allocate_passenger_pathway),
move_type(move_type)
{
    // initialize helper variables
    pathway_indices = new unsigned int[num_driver_pathways];
    for (unsigned int i = 0; i < num_driver_pathways; i++) {
        pathway_indices[i] = i;
    }
    
    for (size_t i = 0; i < 2; i++) {
        move_probs.push_back(0);
        move_log_liks.push_back(0);
    }
    
    if (move_type == GIBBS) {
        size_t num_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
        for (size_t k = 0; k < num_pathways; k++) {
            gibbs_log_liks.push_back(0);
            gibbs_probs.push_back(0);
        }
    }
}

void LinearProgressionModelCache::set_initial_state(gsl_rng *random, LinearProgressionState *prev_state)
{
    this->initial_state = prev_state->increase_num_drivers(random);
}

unsigned long LinearProgressionModelCache::num_iterations()
{
    return num_smc_iter;
}

double LinearProgressionModelCache::get_temperature(size_t t)
{
    return (double)t/num_smc_iter;
}

shared_ptr<LinearProgressionState> LinearProgressionModelCache::propose_initial(gsl_rng *random, double &log_w, LinearProgressionParameters &params)
{
    if (params.get_bgp() != curr_bgp || params.get_fbp() != curr_fbp) {
        reset_cache();
        curr_fbp = params.get_fbp();
        curr_bgp = params.get_bgp();
    }
    
    if (initial_state == 0) {
        if (num_smc_iter == 1) {
            // basically importance sampling -- sample from prior
            LinearProgressionState *state = new LinearProgressionState(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);
            log_w = compute_pathway_likelihood_cache_version2(obs, row_sum, *state, params, _dict_active_probs, _dict_inactive_probs);
            state->set_log_lik(log_w);
            shared_ptr<LinearProgressionState> ret(state);
            return ret;
        }

        LinearProgressionState state(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);
        state.sample_from_prior(random);
        double log_lik = compute_pathway_likelihood_cache_version2(obs, row_sum, state, params, _dict_active_probs, _dict_inactive_probs);
        state.set_log_lik(log_lik);
        return propose_next(random, 0, state, log_w, params);
    } else {
        // sample a state by placing one gene at random from initial state to a new pathway
        // copy state
        double log_lik = compute_pathway_likelihood_cache_version2(obs, row_sum, *initial_state, params, _dict_active_probs, _dict_inactive_probs);
        initial_state->set_log_lik(log_lik);
        return propose_next(random, 0, *initial_state, log_w, params);
    }
}

shared_ptr<LinearProgressionState>LinearProgressionModelCache::propose_next(gsl_rng *random, int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params)
{
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    log_w = temperature_diff * curr.get_log_lik();
    
    // make a copy
    LinearProgressionState *new_state = new LinearProgressionState(curr);
    if (move_type == GIBBS) {
        gibbs_kernel(random, t, *new_state, params);
    } else if (move_type == MH){
        mh_kernel(random, t, *new_state, params);
    } else {
        cerr << "Error: Invalid move type provided." << endl;
        exit(-1);
    }
    shared_ptr<LinearProgressionState> ret(new_state);
    return ret;
}

void LinearProgressionModelCache::swap_pathway_move(gsl_rng *random, LinearProgressionState &state, LinearProgressionParameters &params)
{
    double prev_log_lik = state.get_log_lik();
    move_log_liks[0] = prev_log_lik;
    
    // shuffle indices: swap the first two pathways
    gsl_ran_shuffle(random, pathway_indices, num_driver_pathways, sizeof(unsigned int));
    unsigned int pathway1 = pathway_indices[0];
    unsigned int pathway2 = pathway_indices[1];
    state.swap_pathways(pathway1, pathway2);
    move_log_liks[1] = compute_pathway_likelihood_cache_version2(obs, row_sum, state, params, _dict_active_probs, _dict_inactive_probs);
    normalize(move_log_liks, move_probs);
    size_t k = multinomial(random, move_probs);
    if (k == 0) {
        // revert
        state.swap_pathways(pathway2, pathway1);
    }
    state.set_log_lik(move_log_liks[k]);
}

void LinearProgressionModelCache::mh_kernel(gsl_rng *random, int t, LinearProgressionState &state, LinearProgressionParameters &params)
{
    
    for (size_t i = 0; i < num_mcmc_iter; i++) {
        if (num_driver_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < pathway_swap_prob) {
                swap_pathway_move(random, state, params);
                continue;
            }
        }
        
        // sample a gene at random and sample pathway at random
        double prev_log_lik = state.get_log_lik();
        size_t gene_idx = gsl_rng_uniform_int(random, num_genes);
        size_t old_pathway = state.get_pathway_membership_of(gene_idx);
        size_t new_pathway = gsl_rng_uniform_int(random, state.get_num_pathways());
        state.update_pathway_membership(gene_idx, new_pathway);
        double new_log_lik = compute_pathway_likelihood_cache_version2(obs, row_sum, state, params, _dict_active_probs, _dict_inactive_probs);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - prev_log_lik)) {
            // accept
            state.set_log_lik(new_log_lik);
        } else {
            // revert
            state.update_pathway_membership(gene_idx, old_pathway);
        }
    }
}

void LinearProgressionModelCache::gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &state, LinearProgressionParameters &params)
{
    for (size_t i = 0; i < num_mcmc_iter; i++) {
        if (num_driver_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < pathway_swap_prob) {
                swap_pathway_move(random, state, params);
                continue;
            }
        }
        
        // sample a gene at random and sample pathway at random
        double prev_log_lik = state.get_log_lik();
        size_t gene_idx = gsl_rng_uniform_int(random, num_genes);
        size_t old_pathway = state.get_pathway_membership_of(gene_idx);
        for (size_t k = 0; k < state.get_num_pathways(); k++) {
            if (k == old_pathway) {
                gibbs_log_liks[k] = prev_log_lik;
            }
            state.update_pathway_membership(gene_idx, k);
            gibbs_log_liks[k] = compute_pathway_likelihood_cache_version2(obs, row_sum, state, params, _dict_active_probs, _dict_inactive_probs);
        }
        normalize(gibbs_log_liks, gibbs_probs);
        size_t new_pathway = multinomial(random, gibbs_probs);
        state.update_pathway_membership(gene_idx, new_pathway);
        state.set_log_lik(gibbs_log_liks[new_pathway]);
    }
}

void LinearProgressionModelCache::reset_cache()
{
    _dict_active_probs.clear();
    _dict_inactive_probs.clear();
}

LinearProgressionModelCache::~LinearProgressionModelCache()
{
    delete [] pathway_indices;
    delete initial_state;
}
