//
//  lpm_model.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <assert.h>
#include <cmath>
#include <chrono>
#include <numeric>
#include <iostream>
#include <list>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "lpm_model.hpp"

double log_uniform_prior;

bool perform_swap_move(unsigned int num_driver_pathways, double unif, double pathway_swap_prob)
{
    return (num_driver_pathways >= 2 && unif < pathway_swap_prob);
}

LinearProgressionModel::LinearProgressionModel(unsigned int num_genes,
                                               unsigned int num_driver_pathways,
                                               unsigned int num_smc_iter,
                                               const gsl_matrix &obs,
                                               const vector<unsigned int> &row_sum,
                                               double pathway_swap_prob,
                                               bool allocate_passenger_pathway,
                                               bool is_lik_tempered) :
LinearProgressionModel(num_genes, num_driver_pathways, num_smc_iter, 0, obs, row_sum, pathway_swap_prob, allocate_passenger_pathway, is_lik_tempered)
{
}

LinearProgressionModel::LinearProgressionModel(unsigned int num_genes,
                                               unsigned int num_driver_pathways,
                                               unsigned int num_smc_iter,
                                               unsigned int num_mcmc_iter,
                                               const gsl_matrix &obs,
                                               const vector<unsigned int> &row_sum,
                                               double pathway_swap_prob,
                                               bool allocate_passenger_pathway,
                                               bool is_lik_tempered) :
num_genes(num_genes), num_driver_pathways(num_driver_pathways), num_smc_iter(num_smc_iter),
num_mcmc_iter(num_mcmc_iter), obs(obs), row_sum(row_sum), pathway_swap_prob(pathway_swap_prob),
allocate_passenger_pathway(allocate_passenger_pathway), is_lik_tempered(is_lik_tempered)
{
    if (num_mcmc_iter <= 0) {
        cerr << "Error: number of MCMC kernel iterations must be greater than 0." << endl;
        exit(-1);
    }
    // initialize helper variables
    pathway_indices = new unsigned int[num_driver_pathways];
    for (unsigned int i = 0; i < num_driver_pathways; i++) {
        pathway_indices[i] = i;
    }

    log_uniform_prior = -log_pathway_uniform_prior(num_driver_pathways, num_genes, allocate_passenger_pathway);
}

unsigned long LinearProgressionModel::num_iterations()
{
    return num_smc_iter;
}

double LinearProgressionModel::get_temperature(unsigned int t)
{
    return (double)(t)/num_smc_iter;
}

shared_ptr<LinearProgressionState> LinearProgressionModel::propose_initial(gsl_rng *random, double &log_w, LinearProgressionParameters &params)
{
    double log_lik = DOUBLE_NEG_INF;
    LinearProgressionState *state = 0;
    if (prev_pop == 0) {
        state = new LinearProgressionState(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);
        state->sample_pathway(random);
   } else {
        // sample from the empirical distribution or q_0
        double unif = gsl_ran_flat(random, 0.0, 1.0);
        if (unif < alpha) {
            state = new LinearProgressionState(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);
            state->sample_pathway(random);
        } else {
            // sample from empirical distribution
            double sum = 0.0;
            double w = 0.0;
            unif = gsl_ran_flat(random, 0.0, 1.0);
            for (auto it = prev_pop->begin(); it != prev_pop->end(); ++it) {
                w = (double)it->second/total_mass;
                if (unif < sum + w) {
                    state = new LinearProgressionState(it->first);
                    break;
                }
                sum += w;
            }
            assert(state != 0);
        }
    }
    log_lik = compute_pathway_likelihood(*state, params);
    state->set_log_lik(log_lik);
    log_w = initial_log_weight(*state, params, false);

    shared_ptr<LinearProgressionState> ret(state);
    if (log_lik == DOUBLE_NEG_INF) {
        cerr << "Error: log lik is neg inf! This is a bug as an empty pathway should not have been sampled." << endl;
        exit(-1);
    }

    return ret;
}

shared_ptr<LinearProgressionState>LinearProgressionModel::propose_next(gsl_rng *random, unsigned int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params)
{
    // make a copy of the curr state
    LinearProgressionState *new_state = new LinearProgressionState(curr);    
    log_w = log_weight_helper(t, curr, curr, params, false);
    mh_kernel(random, t, *new_state, params);
    shared_ptr<LinearProgressionState> ret(new_state);
    return ret;
}

void LinearProgressionModel::mh_kernel(gsl_rng *random, double t, LinearProgressionState &state, LinearProgressionParameters &params)
{
    double log_accept_ratio;
    double prev_log_lik = state.get_log_lik();
    double new_log_lik;
    double temperature = get_temperature(t+1);
    
    unsigned int *pathways_to_swap = new unsigned int[2];
    double unif = 0.0;
    bool is_swap_move = false;
    unsigned gene_idx, pathway1, pathway2;

    if (prev_log_lik == DOUBLE_NEG_INF) {
        cerr << "Error: previous log lik is neg inf! This means empty pathway was sampled in the previous step." << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < num_mcmc_iter; i++) {
        unif = gsl_ran_flat(random, 0.0, 1.0);
        is_swap_move = perform_swap_move(num_driver_pathways, unif, pathway_swap_prob);
        if (is_swap_move) {
            gsl_ran_choose(random, pathways_to_swap, 2, pathway_indices, num_driver_pathways, sizeof(unsigned int));
            pathway1 = pathways_to_swap[0];
            pathway2 = pathways_to_swap[1];
            state.swap_pathways(pathway1, pathway2);
        } else {
            // sample a gene at random and sample pathway at random
            gene_idx = gsl_rng_uniform_int(random, num_genes);
            pathway1 = state.get_pathway_membership_of(gene_idx);
            pathway2 = gsl_rng_uniform_int(random, state.get_num_pathways());
            state.update_pathway_membership(gene_idx, pathway2);
        }

        new_log_lik = compute_pathway_likelihood(state, params);

        if (is_lik_tempered)
            log_accept_ratio = temperature * (new_log_lik - prev_log_lik);
        else
            log_accept_ratio = (new_log_lik - prev_log_lik);
        unif = gsl_ran_flat(random, 0.0, 1.0);
        if (log(unif) < log_accept_ratio) {
            // accept
            state.set_log_lik(new_log_lik);
            prev_log_lik = new_log_lik;
        } else {
            // revert
            if (is_swap_move) {
                state.swap_pathways(pathway2, pathway1);
            } else {
                state.update_pathway_membership(gene_idx, pathway1);
            }
            state.set_log_lik(prev_log_lik);
        }
    }
    
    delete [] pathways_to_swap;
}

void LinearProgressionModel::set_particle_population(const vector<shared_ptr<LinearProgressionState> > &particles)
{
    unsigned int n_particles = particles.size();

    if (prev_pop == 0) {
        prev_pop = new unordered_map<LinearProgressionState, unsigned int, hash<LinearProgressionState> >();
    } else {
        prev_pop->clear();
    }
    
    // normalize the probs and copy the contents
    total_mass = n_particles;
    for (size_t i = 0; i < n_particles; i++) {
        const LinearProgressionState &state = *particles[i].get();
        if (prev_pop->count(state) == 0) {
            (*prev_pop)[state] = 0.0;
        }

        (*prev_pop)[state] += 1;
    }
    cout << "Empirical estimate size: " << prev_pop->size() << endl;
}

double LinearProgressionModel::initial_log_weight(const LinearProgressionState &state, const LinearProgressionParameters &params, bool recompute_log_lik)
{
    double log_weight = 0.0;
    double log_proposal = 0.0;
    double log_lik = recompute_log_lik ? compute_pathway_likelihood(state, params) : state.get_log_lik();

    if (is_lik_tempered) {
        double temperature = get_temperature(1);
        log_weight = temperature * log_lik;
    } else {
        log_weight = log_lik;
    }

    if (prev_pop != 0 && prev_pop->count(state) > 0) {
        log_proposal = log(alpha) + log_pathway_proposal(state, num_genes);
        double w = (double)(*prev_pop)[state] / total_mass;
        log_proposal = log_add(log_proposal, log(1 - alpha) + log(w));
    } else {
        log_proposal = log_pathway_proposal(state, num_genes);
    }
    log_weight += log_uniform_prior;
    log_weight -= log_proposal;
    return log_weight;
}

double LinearProgressionModel::log_weight_helper(unsigned int t, const LinearProgressionState &prev_state, const LinearProgressionState &curr_state, const LinearProgressionParameters &params, bool recompute_log_lik)
{
    double log_weight = 0.0;
    double curr_temp = is_lik_tempered ? get_temperature(t+1) : 1.0;
    double prev_temp = is_lik_tempered ? get_temperature(t) : 1.0;
//    double curr_log_lik = recompute_log_lik ? compute_pathway_likelihood(curr_state, params) : curr_state.get_log_lik();
    double prev_log_lik = recompute_log_lik ? compute_pathway_likelihood(prev_state, params) : prev_state.get_log_lik();
    log_weight = (curr_temp - prev_temp) * prev_log_lik;
    return log_weight;
}

double LinearProgressionModel::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<LinearProgressionState> > &genealogy, const LinearProgressionParameters &params)
{
    if (t == 0) {
        const LinearProgressionState &state = *(genealogy->get_state_ptr_at(t).get());
        return initial_log_weight(state, params, true);
    } else {
        const LinearProgressionState &prev_state = *(genealogy->get_state_ptr_at(t-1).get());
        const LinearProgressionState &curr_state = *(genealogy->get_state_ptr_at(t).get());
        return log_weight_helper(t, prev_state, curr_state, params, true);
    }
}

LinearProgressionModel::~LinearProgressionModel()
{
    delete [] pathway_indices;
    delete prev_pop;
}

