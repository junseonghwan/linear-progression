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
                                               size_t num_iter,
                                               size_t num_mh_iter,
                                               gsl_matrix_view *obs,
                                               bool allocate_passenger_pathway)
{
    this->num_genes = num_genes;
    this->num_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
    this->num_iter = num_iter;
    this->num_mh_iter = num_mh_iter;
    this->obs = obs;
    this->allocate_passenger_pathway = allocate_passenger_pathway;
    
    indices = new unsigned int[num_driver_pathways];
    for (unsigned int i = 0; i < num_driver_pathways; i++) {
        indices[i] = i;
    }
}

LinearProgressionModel::LinearProgressionModel(size_t num_genes,
                       size_t num_driver_pathways,
                       size_t num_iter,
                       size_t num_mh_iter,
                       gsl_matrix_view *obs,
                       LinearProgressionState *initial,
                       bool allocate_passenger_pathway) :
LinearProgressionModel(num_genes, num_driver_pathways, num_iter, num_mh_iter, obs, allocate_passenger_pathway)
{
    this->initial = initial;
    this->initial->set_num_pathways(num_driver_pathways);
    this->initial->set_log_lik(DOUBLE_NEG_INF);
}

void LinearProgressionModel::set_initial_population(ParticlePopulation<LinearProgressionState> *pop)
{
    initial_pop = pop;
}

unsigned long LinearProgressionModel::num_iterations()
{
    return num_iter;
}

double LinearProgressionModel::get_temperature(size_t t)
{
    return (double)t/num_iter;
}

pair<LinearProgressionState, double> *LinearProgressionModel::propose_initial(gsl_rng *random, LinearProgressionParameters &params)
{
    if (initial == 0 && initial_pop == 0) {
        // sample from prior
        LinearProgressionState *state = new LinearProgressionState(num_genes, num_pathways, allocate_passenger_pathway);
        state->initialize(random);
        double log_lik = compute_pathway_likelihood(*obs, *state, params);
        state->set_log_lik(log_lik);

        if (num_iter == 1) {
            // perform importance sampling
            return new pair<LinearProgressionState, double>(*state, log_lik);
        }
        return propose_next(random, 0, *state, params);
    } else if (initial_pop == 0) {
        return propose_next(random, 0, *initial, params);
    } else {
        // prior is provided by initial population
        vector<double> probs(initial_pop->get_num_particles());
        normalize(*(initial_pop->get_log_weights()), probs);
        size_t idx = multinomial(random, probs);
        LinearProgressionState &s = initial_pop->get_particles()->at(idx);
        return propose_next(random, 0, s, params);

    }
}

pair<LinearProgressionState, double> *LinearProgressionModel::propose_next(gsl_rng *random, int t, LinearProgressionState curr, LinearProgressionParameters &params)
{
    //double temperature = get_temperature(t+1);
    double temperature_diff = get_temperature(t+1) - get_temperature(t);
    double log_alpha = temperature_diff * curr.get_log_lik();

    for (size_t i = 0; i < num_mh_iter; i++) {
        double prev_log_lik = curr.get_log_lik();
        
        if (num_pathways >= 2) {
            double u = gsl_ran_flat(random, 0.0, 1.0);
            if (u < 0.2) {
                vector<double> log_liks(2);
                log_liks[0] = prev_log_lik;
                
                // perform pathway swap move
                gsl_ran_shuffle(random, indices, num_pathways, sizeof(unsigned int));
                unsigned int pathway1 = indices[0];
                unsigned int pathway2 = indices[1];
                for (size_t g = 0; g < num_genes; g++) {
                    if (curr.get_pathway().at(g) == pathway1) {
                        curr.update_pathway_membership(g, pathway2);
                    } else if (curr.get_pathway().at(g) == pathway2) {
                        curr.update_pathway_membership(g, pathway1);
                    }
                }
                log_liks[1] = compute_pathway_likelihood(*obs, curr, params);
                vector<double> probs(2);
                normalize(log_liks, probs);
                size_t k = multinomial(random, probs);
                if (k == 0) {
                    // revert
                    for (size_t g = 0; g < num_genes; g++) {
                        if (curr.get_pathway().at(g) == pathway2) {
                            curr.update_pathway_membership(g, pathway1);
                        } else if (curr.get_pathway().at(g) == pathway1) {
                            curr.update_pathway_membership(g, pathway2);
                        }
                    }
                }
                
                continue;
            }
        }
        // sample a gene at random and sample pathway
        size_t gene_idx = gsl_rng_uniform_int(random, num_genes);
        size_t old_pathway = curr.get_pathway().at(gene_idx);
        vector<double> log_liks(num_pathways);
        for (size_t new_pathway = 0; new_pathway < num_pathways; new_pathway++) {
            if (old_pathway == new_pathway) {
                log_liks[new_pathway] = prev_log_lik;
                continue;
            }
            curr.update_pathway_membership(gene_idx, new_pathway);
            double new_log_lik = compute_pathway_likelihood(*obs, curr, params);
            log_liks[new_pathway] = new_log_lik;
        }
        vector<double> probs(num_pathways);
        normalize(log_liks, probs);
        double sum = 0.0;
        for (size_t i = 0; i < probs.size(); i++) {
            sum += probs[i];
        }
        if (abs(sum - 1.0) > 1e-6) {
            cerr << "Error in nomralization code!" << endl;
            exit(-1);
        }
        size_t k = multinomial(random, probs);
        curr.update_pathway_membership(gene_idx, k);
        curr.set_log_lik(log_liks[k]);
    }

    return new pair<LinearProgressionState, double>(curr, log_alpha);
}
