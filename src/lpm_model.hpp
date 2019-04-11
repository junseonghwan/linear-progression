//
//  lpm_model.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#ifndef lpm_model_hpp
#define lpm_model_hpp

#include <unordered_map>

#include <gsl/gsl_rng.h>
#include <spf/smc_model.hpp>
#include <spf/particle_population.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_params.hpp"
#include "lpm_state.hpp"

enum MoveType
{
    GIBBS,
    MH
};

class LinearProgressionModel : public ProblemSpecification<LinearProgressionState, LinearProgressionParameters>
{
    unsigned int num_genes;
    unsigned int num_driver_pathways;
    unsigned int num_smc_iter;
    unsigned int num_mcmc_iter;
    const gsl_matrix &obs;
    const vector<unsigned int> &row_sum;
    double pathway_swap_prob;
    bool allocate_passenger_pathway = false;
    MoveType move_type;

    LinearProgressionState *initial_state = 0;

    // helper variables to limit the number of times new vector/array that are allocated
    unsigned int *pathway_indices;

    vector<double> move_log_liks;
    vector<double> move_probs;

    vector<double> gibbs_log_liks;
    vector<double> gibbs_probs;
    
    double get_temperature(unsigned int t);
    void swap_pathway_move(gsl_rng *random, LinearProgressionState &state, LinearProgressionParameters &params);
    void mh_kernel(gsl_rng *random, LinearProgressionState &curr, LinearProgressionParameters &params);
    void gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &curr, LinearProgressionParameters &params);
    
public:

    LinearProgressionModel(unsigned int num_genes,
                           unsigned int num_driver_pathways,
                           unsigned int num_iter,
                           unsigned int num_mcmc_iter,
                           const gsl_matrix &obs,
                           const vector<unsigned int> &row_sum,
                           double pathway_swap_prob,
                           bool allocate_passenger_pathway,
                           MoveType move_type);
    void set_initial_state(gsl_rng *random, LinearProgressionState *prev_state);
    unsigned long num_iterations();
    shared_ptr<LinearProgressionState> propose_initial(gsl_rng *random, double &log_w, LinearProgressionParameters &params);
    shared_ptr<LinearProgressionState> propose_next(gsl_rng *random, unsigned int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params);
    double log_weight(unsigned int t, const LinearProgressionState &curr, const LinearProgressionParameters &params);
    ~LinearProgressionModel();
    
};

#endif /* lpm_model_hpp */
