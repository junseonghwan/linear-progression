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
    bool is_lik_tempered = false;

    // variables related to initial proposal: (1-alpha) \hat{p}(x|y) + alpha p(x)
    double alpha = 0.5; // mixture proportion for the prior proposal within cSMC context
    // prev_pop to be used for initial proposal
    unordered_map<LinearProgressionState, unsigned int, hash<LinearProgressionState> > *prev_pop = 0;
    unsigned int total_mass = 0;

    // helper variables to limit the number of times new vector/array that are allocated
    unsigned int *pathway_indices;

    double get_temperature(unsigned int t);
    void rw_move(gsl_rng *random, double t, LinearProgressionState &curr, LinearProgressionParameters &params);
    void mh_kernel(gsl_rng *random, double t, LinearProgressionState &curr, LinearProgressionParameters &params);
    double initial_log_weight(const LinearProgressionState &state, const LinearProgressionParameters &params, bool recompute_log_lik);
    double log_weight_helper(unsigned int t, const LinearProgressionState &prev_state, const LinearProgressionState &curr_state, const LinearProgressionParameters &params, bool recompute_log_lik);
public:

    // constructor for random walk proposal
    LinearProgressionModel(unsigned int num_genes,
                           unsigned int num_driver_pathways,
                           unsigned int num_smc_iter,
                           const gsl_matrix &obs,
                           const vector<unsigned int> &row_sum,
                           double pathway_swap_prob,
                           bool allocate_passenger_pathway,
                           bool is_lik_tempered);

    // constructor for MH kernel proposal
    LinearProgressionModel(unsigned int num_genes,
                           unsigned int num_driver_pathways,
                           unsigned int num_smc_iter,
                           unsigned int num_mcmc_iter,
                           const gsl_matrix &obs,
                           const vector<unsigned int> &row_sum,
                           double pathway_swap_prob,
                           bool allocate_passenger_pathway,
                           bool is_lik_tempered);
    
    unsigned long num_iterations() override;
    shared_ptr<LinearProgressionState> propose_initial(gsl_rng *random, double &log_w, LinearProgressionParameters &params) override;
    shared_ptr<LinearProgressionState> propose_next(gsl_rng *random, unsigned int t, const LinearProgressionState &curr, double &log_w, LinearProgressionParameters &params) override;
    double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<LinearProgressionState> > &genealogy, const LinearProgressionParameters &params) override;
    void set_particle_population(const vector<shared_ptr<LinearProgressionState> > &particles) override;
    ~LinearProgressionModel();
    
};

#endif /* lpm_model_hpp */
