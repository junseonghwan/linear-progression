//
//  lpm_likelihood.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#ifndef lpm_likelihood_hpp
#define lpm_likelihood_hpp

#include <gsl/gsl_blas.h>

#include "lpm_params.hpp"
#include "lpm_state.hpp"

// matrix product between m-th row of obs and mutation matrix induced by pathway_membership
void fast_matrix_product(gsl_matrix_view &obs, size_t m, vector<size_t> &pathway_membership, size_t num_pathways, vector<size_t> &ret);

// compute log likelihood for a given number of mutations given pathway size and error probabilities
double compute_log_lik_active(double n_mutations, double pathway_size, double bgp, double fbp);
double compute_log_lik_not_active(double n_mutations, double pathway_size, double bgp, double fbp);
double compute_log_lik_passenger(double n_mutations, double bgp);

double compute_likelihood_for_sample(vector<size_t> &r, vector<size_t> &pathway_sizes, size_t stage, double bgp, double fbp);
double compute_pathway_likelihood(gsl_matrix_view &obs, LinearProgressionState &state, LinearProgressionParameters &params);

//double compute_pathway_likelihood(gsl_matrix_view &Y, gsl_matrix_view &X, vector<size_t> &stages, double fbp, double bgp);

#endif /* lpm_likelihood_hpp */
