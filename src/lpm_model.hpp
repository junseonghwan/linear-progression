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
    size_t num_driver_pathways;
    size_t num_smc_iter;
    size_t num_mcmc_iter;
    gsl_matrix_view &obs;
    vector<size_t> &row_sum;
    double pathway_swap_prob;
    bool allocate_passenger_pathway = false;
    bool use_gibbs_kernel = false;
    
    //ParticlePopulation<LinearProgressionState> *initial_pop = 0;
    LinearProgressionState *initial_state = 0;
    
    // helper variables to limit the number of times new vector/array are allocated
    unsigned int *pathway_indices;

    vector<double> move_log_liks;
    vector<double> move_probs;
    
    vector<double> gibbs_log_liks;
    vector<double> gibbs_probs;
    
    double get_temperature(size_t t);
    pair<LinearProgressionState, double> *mh_kernel(gsl_rng *random, int t, LinearProgressionState &curr,  LinearProgressionParameters &params);
    pair<LinearProgressionState, double> *gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &curr,  LinearProgressionParameters &params);


public:
    LinearProgressionModel(size_t num_genes,
                           size_t num_driver_pathways,
                           size_t num_iter,
                           size_t num_mcmc_iter,
                           gsl_matrix_view &obs,
                           vector<size_t> &row_sum,
                           double pathway_swap_prob = 0.2,
                           bool allocate_passenger_pathway = true,
                           bool use_gibbs_kernel = true);
    //void set_initial_population(ParticlePopulation<LinearProgressionState> *pop);
    void set_initial_state(LinearProgressionState *prev_state);
    unsigned long num_iterations();
    pair<LinearProgressionState, double> *propose_initial(gsl_rng *random, LinearProgressionParameters &params);
    pair<LinearProgressionState, double> *propose_next(gsl_rng *random, int t, LinearProgressionState curr, LinearProgressionParameters &params);
    ~LinearProgressionModel() { }
};

#endif /* lpm_model_hpp */
