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

#include "lpm_model.hpp"

LinearProgressionModel::LinearProgressionModel(size_t num_genes,
                                               size_t num_driver_pathways,
                                               size_t num_smc_iter,
                                               size_t num_mcmc_iter,
                                               gsl_matrix_view &obs,
                                               vector<size_t> &row_sum,
                                               double pathway_swap_prob,
                                               bool allocate_passenger_pathway,
                                               bool use_gibbs_kernel) :
num_genes(num_genes), num_driver_pathways(num_driver_pathways), num_smc_iter(num_smc_iter), num_mcmc_iter(num_mcmc_iter),
obs(obs), row_sum(row_sum), pathway_swap_prob(pathway_swap_prob), allocate_passenger_pathway(allocate_passenger_pathway),
use_gibbs_kernel(use_gibbs_kernel)
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
    
    if (use_gibbs_kernel) {
        size_t num_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
        for (size_t k = 0; k < num_pathways; k++) {
            gibbs_log_liks.push_back(0);
            gibbs_probs.push_back(0);
        }
    }
}

//void LinearProgressionModel::set_initial_population(ParticlePopulation<LinearProgressionState> *pop)
//{
//    initial_pop = pop;
//}

void LinearProgressionModel::set_initial_state(LinearProgressionState *prev_state)
{
    this->initial_state = prev_state->increase_num_drivers();
}

unsigned long LinearProgressionModel::num_iterations()
{
    return num_smc_iter;
}

double LinearProgressionModel::get_temperature(size_t t)
{
    return (double)t/num_smc_iter;
}

pair<LinearProgressionState, double> *LinearProgressionModel::propose_initial(gsl_rng *random, LinearProgressionParameters &params)
{
    if (initial_state == 0) {
        // sample from prior
        LinearProgressionState *state = new LinearProgressionState(num_genes, num_driver_pathways, allocate_passenger_pathway);
        //state->sample_min_valid_pathway(random);
        state->sample_from_prior(random);
        double log_lik = compute_pathway_likelihood(obs, row_sum, *state, params);
        state->set_log_lik(log_lik);

        if (num_smc_iter == 1) {
            // perform importance sampling
            return new pair<LinearProgressionState, double>(*state, log_lik);
        }
        return propose_next(random, 0, *state, params);
    } else {
        // sample a state by placing one gene at random from initial state to a new pathway
        // copy state
        double log_lik = compute_pathway_likelihood(obs, row_sum, *initial_state, params);
        initial_state->set_log_lik(log_lik);
        return propose_next(random, 0, *initial_state, params);
    }
}

pair<LinearProgressionState, double> *LinearProgressionModel::propose_next(gsl_rng *random, int t, LinearProgressionState curr, LinearProgressionParameters &params)
{
    if (use_gibbs_kernel) {
        return gibbs_kernel(random, t, curr, params);
    } else {
        return mh_kernel(random, t, curr, params);
    }
}

pair<LinearProgressionState, double> *LinearProgressionModel::mh_kernel(gsl_rng *random, int t, LinearProgressionState &curr,  LinearProgressionParameters &params)
{
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    double log_alpha = temperature_diff * curr.get_log_lik();
    
    for (size_t i = 0; i < num_mcmc_iter; i++) {
        double prev_log_lik = curr.get_log_lik();
        
        if (num_driver_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < pathway_swap_prob) {
                move_log_liks[0] = prev_log_lik;
                
                // shuffle indices: swap the first two pathways
                gsl_ran_shuffle(random, pathway_indices, num_driver_pathways, sizeof(unsigned int));
                unsigned int pathway1 = pathway_indices[0];
                unsigned int pathway2 = pathway_indices[1];
                curr.swap_pathways(pathway1, pathway2);
                move_log_liks[1] = compute_pathway_likelihood(obs, row_sum, curr, params);
                normalize(move_log_liks, move_probs);
                size_t k = multinomial(random, move_probs);
                if (k == 0) {
                    // revert
                    curr.swap_pathways(pathway2, pathway1);
                }
                curr.set_log_lik(move_log_liks[k]);
                continue;
            }
        }
        
        // sample a gene at random and sample pathway at random
        size_t gene_idx = gsl_rng_uniform_int(random, num_genes);
        size_t old_pathway = curr.get_pathway_membership(gene_idx);
        size_t new_pathway = gsl_rng_uniform_int(random, curr.get_num_pathways());
        curr.update_pathway_membership(gene_idx, new_pathway);
        double new_log_lik = compute_pathway_likelihood(obs, row_sum, curr, params);
        double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
        if (log_u < (new_log_lik - prev_log_lik)) {
            // accept
            curr.set_log_lik(new_log_lik);
        } else {
            // revert
            curr.update_pathway_membership(gene_idx, old_pathway);
        }
    }
    return new pair<LinearProgressionState, double>(curr, log_alpha);
}

pair<LinearProgressionState, double> *LinearProgressionModel::gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &curr,  LinearProgressionParameters &params)
{
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    double log_alpha = temperature_diff * curr.get_log_lik();
    
    for (size_t i = 0; i < num_mcmc_iter; i++) {
        double prev_log_lik = curr.get_log_lik();
        
        if (num_driver_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < pathway_swap_prob) {
                move_log_liks[0] = prev_log_lik;
                
                // shuffle indices: swap the first two pathways
                gsl_ran_shuffle(random, pathway_indices, num_driver_pathways, sizeof(unsigned int));
                unsigned int pathway1 = pathway_indices[0];
                unsigned int pathway2 = pathway_indices[1];
                curr.swap_pathways(pathway1, pathway2);
                move_log_liks[1] = compute_pathway_likelihood(obs, row_sum, curr, params);
                normalize(move_log_liks, move_probs);
                size_t k = multinomial(random, move_probs);
                if (k == 0) {
                    // revert
                    curr.swap_pathways(pathway2, pathway1);
                }
                curr.set_log_lik(move_log_liks[k]);
                continue;
            }
        }

        // sample a gene at random and sample pathway at random
        size_t gene_idx = gsl_rng_uniform_int(random, num_genes);
        size_t old_pathway = curr.get_pathway_membership(gene_idx);
        for (size_t k = 0; k < curr.get_num_pathways(); k++) {
            if (k == old_pathway) {
                gibbs_log_liks[k] = prev_log_lik;
            }
            curr.update_pathway_membership(gene_idx, k);
            gibbs_log_liks[k] = compute_pathway_likelihood(obs, row_sum, curr, params);
        }
        normalize(gibbs_log_liks, gibbs_probs);
        size_t new_pathway = multinomial(random, gibbs_probs);
        curr.update_pathway_membership(gene_idx, new_pathway);
        curr.set_log_lik(gibbs_log_liks[new_pathway]);
    }
    return new pair<LinearProgressionState, double>(curr, log_alpha);
}
