//
//  lpm_model.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <assert.h>
#include <math.h>
#include <numeric>
#include <iostream>
#include <list>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "lpm_model.hpp"

bool perform_swap_move(unsigned int num_driver_pathways, double unif, double pathway_swap_prob)
{
    return (num_driver_pathways >= 2 && unif < pathway_swap_prob);
}

LinearProgressionModel::LinearProgressionModel(unsigned int num_genes,
                                               unsigned int num_driver_pathways,
                                               unsigned int num_smc_iter,
                                               unsigned int num_mcmc_iter,
                                               const gsl_matrix &obs,
                                               const vector<unsigned int> &row_sum,
                                               double pathway_swap_prob,
                                               bool allocate_passenger_pathway) :
num_genes(num_genes), num_driver_pathways(num_driver_pathways), num_smc_iter(num_smc_iter),
num_mcmc_iter(num_mcmc_iter), obs(obs), row_sum(row_sum), pathway_swap_prob(pathway_swap_prob),
allocate_passenger_pathway(allocate_passenger_pathway)
{
    // initialize helper variables
    pathway_indices = new unsigned int[num_driver_pathways];
    for (unsigned int i = 0; i < num_driver_pathways; i++) {
        pathway_indices[i] = i;
    }
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
    double log_lik;
    LinearProgressionState *state = 0;

    if (prev_pop == 0) {

        state = new LinearProgressionState(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);

        state->sample_from_prior(random);
        log_lik = compute_pathway_likelihood(*state, params);
        state->set_log_lik(log_lik);
        shared_ptr<LinearProgressionState> ret(state);

        if (num_smc_iter == 1) {
            log_w = log_lik;
        } else {
            log_w = get_temperature(1) * log_lik;
        }
        return ret;
    }

    // proposal from (1-\alpha) \hat{p}(x|y) + \alpha p(x)
    double unif = gsl_ran_flat(random, 0.0, 1.0);
    if (unif < alpha) {
        // sample from the prior
        state = new LinearProgressionState(obs, row_sum, num_genes, num_driver_pathways, allocate_passenger_pathway);
        state->sample_from_prior(random);
    } else {
        // sample from prev_pop
        unif = gsl_ran_flat(random, 0.0, 1.0);
        double sum = 0.0;
        for (auto it = prev_pop->begin(); it != prev_pop->end(); ++it) {
            if (unif < (sum + it->second)) {
                state = new LinearProgressionState(it->first);
                break;
            }
            sum += it->second;
        }
    }
    assert(state != 0);
    log_lik = compute_pathway_likelihood(*state, params);
    state->set_log_lik(log_lik);
    log_w = log_weight_helper(0, *state, params, false);
    shared_ptr<LinearProgressionState> ret(state);
    return ret;
}

shared_ptr<LinearProgressionState>LinearProgressionModel::propose_next(gsl_rng *random, unsigned int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params)
{
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    log_w = temperature_diff * curr.get_log_lik();

    // make a copy of the curr state
    LinearProgressionState *new_state = new LinearProgressionState(curr);
    mh_kernel(random, t, *new_state, params);
    shared_ptr<LinearProgressionState> ret(new_state);
    return ret;
}

void LinearProgressionModel::mh_kernel(gsl_rng *random, double t, LinearProgressionState &state, LinearProgressionParameters &params)
{
    double log_accept_ratio;
    double prev_log_lik = state.get_log_lik();
    double prev_log_prior = log_pathway_prior(state, num_genes);
    double new_log_lik;
    double new_log_prior;
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
            gsl_ran_sample(random, pathways_to_swap, 2, pathway_indices, num_driver_pathways, sizeof(unsigned int));
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
        new_log_prior = log_pathway_prior(state, num_genes);

        log_accept_ratio = temperature * (new_log_lik - prev_log_lik);
        log_accept_ratio += (new_log_prior - prev_log_prior);
        unif = gsl_ran_flat(random, 0.0, 1.0);
        if (log(unif) < log_accept_ratio) {
            // accept
            state.set_log_lik(new_log_lik);
            prev_log_lik = new_log_lik;
            prev_log_prior = new_log_prior;
        } else {
            // revert
            if (is_swap_move) {
                state.swap_pathways(pathway2, pathway1);
            } else {
                state.update_pathway_membership(gene_idx, pathway1);
            }
        }
    }
}

void LinearProgressionModel::set_particle_population(const vector<vector<shared_ptr<LinearProgressionState> > > &particles, const vector<vector<double> > &log_weights, const vector<double> &log_norms)
{
    unsigned int R = num_iterations();
    // take the samples from the last iteration to update the initial proposal distribution
    const vector<shared_ptr<LinearProgressionState> > &particles_at_termination = particles.at(R - 1);
    const vector<double> &log_weights_at_termination = log_weights.at(R-1);
    double log_norm = log_norms[R-1];

    unsigned int n_particles = particles_at_termination.size();

    if (prev_pop == 0) {
        prev_pop = new unordered_map<LinearProgressionState, double, hash<LinearProgressionState> >();
    } else {
        prev_pop->clear();
    }

    // normalize the probs and copy the contents
    double sum = 0.0;
    double prob;
    for (size_t i = 0; i < n_particles; i++) {
        const LinearProgressionState &state = *particles_at_termination[i].get();
        if (prev_pop->count(state) == 0) {
            (*prev_pop)[state] = 0.0;
        }
        
        prob = exp(log_weights_at_termination[i] - log_norm);
        (*prev_pop)[state] += prob;
        sum += prob;
    }

    assert(abs(sum - 1.0) < 1e-6);
}

double LinearProgressionModel::log_weight_helper(unsigned int t, const LinearProgressionState &state, const LinearProgressionParameters &params, bool recompute_log_lik)
{
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    double log_lik = recompute_log_lik ? compute_pathway_likelihood(state, params) : state.get_log_lik();
    double log_weight = temperature_diff * log_lik;
    if (t == 0) {
        double log_prior = log_pathway_prior(state, num_genes);
        double log_proposal = log(alpha) + log_prior;
        if (prev_pop->count(state) > 0) {
            log_proposal = log_add(log_proposal, log(1 - alpha) + log((*prev_pop)[state]));
        }
        log_weight += log_prior;
        log_weight -= log_proposal;
    }
    return log_weight;
}

double LinearProgressionModel::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<LinearProgressionState> > &genealogy, const LinearProgressionParameters &params)
{
    if (t == 0) {
        const LinearProgressionState &state = *(genealogy->get_state_ptr_at(t).get());
        return log_weight_helper(t, state, params, true);
    } else {
        const LinearProgressionState &prev_state = *(genealogy->get_state_ptr_at(t-1).get());
        return log_weight_helper(t, prev_state, params, true);
    }
}

LinearProgressionModel::~LinearProgressionModel()
{
    delete [] pathway_indices;
    delete prev_pop;
}

//void LinearProgressionModel::swap_pathway_move(gsl_rng *random, double t, LinearProgressionState &state, LinearProgressionParameters &params)
//{
//    double prev_log_lik = state.get_log_lik();
//    double prev_log_prior = log_pathway_prior(state, num_genes);
//    double log_acceptance_ratio;
//
//    // shuffle indices: swap the first two pathways
//    gsl_ran_shuffle(random, pathway_indices, num_driver_pathways, sizeof(unsigned int));
//    unsigned int pathway1 = pathway_indices[0];
//    unsigned int pathway2 = pathway_indices[1];
//    state.swap_pathways(pathway1, pathway2);
//    double new_log_lik = compute_pathway_likelihood(state, params);
//    double new_log_prior = log_pathway_prior(state, num_genes);
//
//    // perform accept-reject
//    double log_unif = log(gsl_ran_flat(random, 0.0, 1.0));
//    log_acceptance_ratio = temperature * (new_log_lik - prev_log_lik);
//    log_acceptance_ratio += (new_log_prior - prev_log_prior);
//    if (log_unif < log_acceptance_ratio) {
//        state.set_log_lik(new_log_lik);
//    } else {
//        // revert
//        state.swap_pathways(pathway2, pathway1);
//    }
//}
//
//void LinearProgressionModel::gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &state, LinearProgressionParameters &params)
//{
//    double temperature = get_temperature(t+1);
//    for (unsigned int i = 0; i < num_mcmc_iter; i++) {
//        if (num_driver_pathways >= 2) {
//            double u = gsl_ran_flat(random, 0.0, 1.0);
//            if (u < pathway_swap_prob) {
//                swap_pathway_move(random, temperature, state, params);
//                continue;
//            }
//        }
//
//        // sample a gene at random and sample pathway at random
//        double prev_log_lik = state.get_log_lik();
//        unsigned int gene_idx = gsl_rng_uniform_int(random, num_genes);
//        unsigned int old_pathway = state.get_pathway_membership_of(gene_idx);
//        for (unsigned int k = 0; k < state.get_num_pathways(); k++) {
//            if (k == old_pathway) {
//                gibbs_log_liks[k] = prev_log_lik;
//            }
//            state.update_pathway_membership(gene_idx, k);
//            gibbs_log_liks[k] = compute_pathway_likelihood(state, params);
//        }
//        normalize(gibbs_log_liks, gibbs_probs);
//        unsigned int new_pathway = multinomial(random, gibbs_probs);
//        state.update_pathway_membership(gene_idx, new_pathway);
//        state.set_log_lik(gibbs_log_liks[new_pathway]);
//    }
//}
