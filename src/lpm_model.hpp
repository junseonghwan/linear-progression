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
    size_t num_genes;
    size_t num_driver_pathways;
    size_t num_smc_iter;
    size_t num_mcmc_iter;
    gsl_matrix &obs;
    vector<size_t> &row_sum;
    double pathway_swap_prob;
    bool allocate_passenger_pathway = false;
    MoveType move_type;

    LinearProgressionState *initial_state = 0;

    // helper variables to limit the number of times new vector/array are allocated
    unsigned int *pathway_indices;

    vector<double> move_log_liks;
    vector<double> move_probs;

    vector<double> gibbs_log_liks;
    vector<double> gibbs_probs;
    
    // cache likelihood calculation
    double cached_bgp, cached_fbp;
    unordered_map<size_t, unordered_map<size_t, double> > _dict_active_probs;
    unordered_map<size_t, unordered_map<size_t, double> > _dict_inactive_probs;
    void initialize_likelihood_table(LinearProgressionParameters &params);
    bool update_cache(LinearProgressionParameters &params);
    //void construct_likelihood_tables(LinearProgressionParameters &params);
    
    double get_temperature(size_t t);
    void swap_pathway_move(gsl_rng *random, LinearProgressionState &state, LinearProgressionParameters &params);
    void mh_kernel(gsl_rng *random, LinearProgressionState &curr, LinearProgressionParameters &params);
    void gibbs_kernel(gsl_rng *random, int t, LinearProgressionState &curr, LinearProgressionParameters &params);
    
public:

    LinearProgressionModel(size_t num_genes,
                           size_t num_driver_pathways,
                           size_t num_iter,
                           size_t num_mcmc_iter,
                           gsl_matrix &obs,
                           vector<size_t> &row_sum,
                           double pathway_swap_prob,
                           bool allocate_passenger_pathway,
                           MoveType move_type);
    void set_initial_state(gsl_rng *random, LinearProgressionState *prev_state);
    unsigned long num_iterations();
    shared_ptr<LinearProgressionState> propose_initial(gsl_rng *random, double &log_w, LinearProgressionParameters &params);
    shared_ptr<LinearProgressionState> propose_next(gsl_rng *random, int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params);
    ~LinearProgressionModel();
    
};

#endif /* lpm_model_hpp */
