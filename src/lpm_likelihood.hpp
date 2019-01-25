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

double compute_likelihood(gsl_matrix_view &ret, size_t m, size_t stage, vector<size_t> pathway_sizes, double bgp, double fbp);
double compute_log_lik(double n_mutations, double pathway_size, double bgp, double fbp);
double compute_log_lik_bg(double n_mutations, double pathway_size, double bgp, double fbp);
double compute_pathway_likelihood(gsl_matrix_view &Y, gsl_matrix_view &X, vector<size_t> &stages, double fbp, double bgp);
double compute_pathway_likelihood(gsl_matrix_view &obs, LinearProgressionState &state, LinearProgressionParameters &params);

#endif /* lpm_likelihood_hpp */
