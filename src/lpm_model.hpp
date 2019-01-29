//
//  lpm_model.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#ifndef lpm_model_hpp
#define lpm_model_hpp

#include <gsl/gsl_rng.h>
#include <spf/smc_model.hpp>
#include <spf/particle_population.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_params.hpp"
#include "lpm_state.hpp"

class LinearProgressionModel : public ProblemSpecification<LinearProgressionState, LinearProgressionParameters>
{
    size_t num_genes;
    size_t num_pathways;
    size_t num_iter;
    size_t num_mcmc_iter;
    gsl_matrix_view *obs = 0;
    double pathway_swap_prob;
    ParticlePopulation<LinearProgressionState> *initial_pop = 0;
    bool allocate_passenger_pathway;
    double get_temperature(size_t t);
    
    unsigned int *indices;

    vector<double> gibbs_log_liks;
    vector<double> gibbs_probs;
    vector<double> swap_log_liks;
    vector<double> swap_probs;

public:
    LinearProgressionModel(size_t num_genes,
                           size_t num_driver_pathways,
                           size_t num_iter,
                           size_t num_mcmc_iter,
                           gsl_matrix_view *obs,
                           double pathway_swap_prob = 0.2,
                           bool allocate_passenger_pathway = false);
    void set_initial_population(ParticlePopulation<LinearProgressionState> *pop);
    unsigned long num_iterations();
    pair<LinearProgressionState, double> *propose_initial(gsl_rng *random, LinearProgressionParameters &params);
    pair<LinearProgressionState, double> *propose_next(gsl_rng *random, int t, LinearProgressionState curr, LinearProgressionParameters &params);
    ~LinearProgressionModel() { }
};

#endif /* lpm_model_hpp */
