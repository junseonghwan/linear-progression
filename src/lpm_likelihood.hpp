//
//  lpm_likelihood.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#ifndef lpm_likelihood_hpp
#define lpm_likelihood_hpp

#include <unordered_map>

#include <gsl/gsl_blas.h>

#include "lpm_params.hpp"
#include "lpm_state.hpp"

// matrix product between m-th row of obs and mutation matrix induced by pathway_membership
void fast_matrix_product(gsl_matrix &obs, unsigned int m, const vector<unsigned int> &pathway_membership, vector<unsigned int> &ret);

double log_pathway_uniform_prior(unsigned int n_pathways, unsigned int n_genes);
double log_pathway_proposal(const LinearProgressionState &pathway, unsigned int n_genes);
double log_pathway_proposal(const vector<unsigned int> &pathway, unsigned int n_pathways);
double log_min_valid_pathway_proposal(const LinearProgressionState &pathway, unsigned int n_genes);

// compute log likelihood for a given number of mutations given pathway size and error probabilities
double compute_log_lik_active(unsigned int n_mutations, unsigned int pathway_size, double bgp, double fbp);
double compute_log_lik_inactive(unsigned int n_mutations, unsigned int pathway_size, double bgp);
//double compute_likelihood_for_sample(const vector<unsigned short> &r,
//                                     vector<unsigned int> &pathway_sizes,
//                                     unsigned int stage,
//                                     double bgp, double fbp);
double compute_likelihood_for_sample(const vector<unsigned short> &r,
                                     const LinearProgressionState &state,
                                     unsigned int stage,
                                     double bgp, double fbp);
double compute_pathway_likelihood(const LinearProgressionState &state,
                                   const LinearProgressionParameters &params);


#endif /* lpm_likelihood_hpp */
