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

LinearProgressionModel::LinearProgressionModel(unsigned int num_genes,
                                               unsigned int num_driver_pathways,
                                               unsigned int num_smc_iter,
                                               unsigned int num_mcmc_iter,
                                               const gsl_matrix &obs,
                                               const vector<unsigned int> &row_sum,
                                               double pathway_swap_prob,
                                               bool allocate_passenger_pathway,
                                               MoveType move_type) :
num_genes(num_genes), num_driver_pathways(num_driver_pathways), num_smc_iter(num_smc_iter),
num_mcmc_iter(num_mcmc_iter), obs(obs), row_sum(row_sum), pathway_swap_prob(pathway_swap_prob),
allocate_passenger_pathway(allocate_passenger_pathway),
move_type(move_type)
{
    // initialize helper variables
    pathway_indices = new unsigned int[num_driver_pathways];
    for (unsigned int i = 0; i < num_driver_pathways; i++) {
        pathway_indices[i] = i;
    }
    
    for (unsigned int i = 0; i < 2; i++) {
        move_probs.push_back(0);
        move_log_liks.push_back(0);
    }

    if (move_type == GIBBS) {
        unsigned int num_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
        for (unsigned int k = 0; k < num_pathways; k++) {
            gibbs_log_liks.push_back(0);
            gibbs_probs.push_back(0);
        }
    }
}

void LinearProgressionModel::set_initial_state(gsl_rng *random, LinearProgressionState *prev_state)
{
    // sample a state by placing one gene at random from initial state to a new pathway
    this->initial_state = prev_state->increase_num_drivers(random);
}

unsigned long LinearProgressionModel::num_iterations()
{
    return num_smc_iter;
}

double LinearProgressionModel::get_temperature(unsigned int t)
{
    return (double)(t)/(num_smc_iter-1);
}

shared_ptr<LinearProgressionState> LinearProgressionModel::propose_initial(gsl_rng *random, double &log_w, LinearProgressionParameters &params)
{
    if (initial_state == 0) {
        LinearProgressionState *state = new LinearProgressionState(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);

        if (!allocate_passenger_pathway) {
            state->sample_from_prior(random);
        } else {
            // in this case, we randomly allocate one gene to driver pathway and the rest to passengers
            state->sample_min_valid_pathway(random);
        }

        double log_lik = compute_pathway_likelihood(*state, params);
        state->set_log_lik(log_lik);
        shared_ptr<LinearProgressionState> ret(state);

        if (num_smc_iter == 1) {
            log_w = log_lik;
        } else {
            log_w = 0.0;
        }

        if (allocate_passenger_pathway) {
            // account for the difference between the prior and the proposal
            // this factor is same for all state so leave it for now:
            // implement this for testing purposes.
        }
        return ret;

    } else {
        double log_lik = compute_pathway_likelihood(*initial_state, params);
        initial_state->set_log_lik(log_lik);
        return propose_next(random, 0, *initial_state, log_w, params);
    }
}

shared_ptr<LinearProgressionState>LinearProgressionModel::propose_next(gsl_rng *random, unsigned int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params)
{
    double temperature_diff = get_temperature(t) - get_temperature(t-1);
    log_w = temperature_diff * curr.get_log_lik();

    // make a copy
    LinearProgressionState *new_state = new LinearProgressionState(curr);
    if (move_type == GIBBS) {
        gibbs_kernel(random, t, *new_state, params);
    } else if (move_type == MH){
        mh_kernel(random, *new_state, params);
    } else {
        cerr << "Error: Invalid move type provided." << endl;
        exit(-1);
    }
    shared_ptr<LinearProgressionState> ret(new_state);
    return ret;
}

double LinearProgressionModel::log_weight(unsigned int t, const LinearProgressionState &state, const LinearProgressionParameters &params)
{
    double temperature_diff = get_temperature(t) - get_temperature(t-1);
    return temperature_diff * compute_pathway_likelihood(state, params);
}

void LinearProgressionModel::swap_pathway_move(gsl_rng *random, LinearProgressionState &state, LinearProgressionParameters &params)
{
    double prev_log_lik = state.get_log_lik();
    move_log_liks[0] = prev_log_lik;
    if (prev_log_lik == DOUBLE_NEG_INF) {
        cout << "neg inf" << endl;
    }
    
    // shuffle indices: swap the first two pathways
    gsl_ran_shuffle(random, pathway_indices, num_driver_pathways, sizeof(unsigned int));
    unsigned int pathway1 = pathway_indices[0];
    unsigned int pathway2 = pathway_indices[1];
    state.swap_pathways(pathway1, pathway2);
    move_log_liks[1] = compute_pathway_likelihood(state, params);
    vector<double> probs(2);
    normalize(move_log_liks, probs);
    unsigned int idx = multinomial(random, probs);
    if (idx == 1) {
        state.set_log_lik(move_log_liks[1]);
    } else {
        // revert
        state.swap_pathways(pathway2, pathway1);
    }
}

void LinearProgressionModel::mh_kernel(gsl_rng *random, LinearProgressionState &state, LinearProgressionParameters &params)
{

    for (unsigned int i = 0; i < num_mcmc_iter; i++) {
        if (num_driver_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < pathway_swap_prob) {
                swap_pathway_move(random, state, params);
                continue;
            }
        }
        
        // sample a gene at random and sample pathway at random
        double prev_log_lik = state.get_log_lik();
        unsigned int gene_idx = gsl_rng_uniform_int(random, num_genes);
        unsigned int old_pathway = state.get_pathway_membership_of(gene_idx);
        unsigned int new_pathway = gsl_rng_uniform_int(random, state.get_num_pathways());
        if (old_pathway == new_pathway) {
            continue;
        }
        state.update_pathway_membership(gene_idx, new_pathway);
        double new_log_lik = compute_pathway_likelihood(state, params);
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

void LinearProgressionModel::gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &state, LinearProgressionParameters &params)
{
    for (unsigned int i = 0; i < num_mcmc_iter; i++) {
        if (num_driver_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < pathway_swap_prob) {
                swap_pathway_move(random, state, params);
                continue;
            }
        }

        // sample a gene at random and sample pathway at random
        double prev_log_lik = state.get_log_lik();
        unsigned int gene_idx = gsl_rng_uniform_int(random, num_genes);
        unsigned int old_pathway = state.get_pathway_membership_of(gene_idx);
        for (unsigned int k = 0; k < state.get_num_pathways(); k++) {
            if (k == old_pathway) {
                gibbs_log_liks[k] = prev_log_lik;
            }
            state.update_pathway_membership(gene_idx, k);
            gibbs_log_liks[k] = compute_pathway_likelihood(state, params);
        }
        normalize(gibbs_log_liks, gibbs_probs);
        unsigned int new_pathway = multinomial(random, gibbs_probs);
        state.update_pathway_membership(gene_idx, new_pathway);
        state.set_log_lik(gibbs_log_liks[new_pathway]);
    }
}

LinearProgressionModel::~LinearProgressionModel()
{
    delete [] pathway_indices;
    delete initial_state;
}
