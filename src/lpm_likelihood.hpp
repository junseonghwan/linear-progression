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

double log_pathway_uniform_prior(unsigned int n_pathways, unsigned int n_genes);
double log_pathway_prior(const LinearProgressionState &pathway, unsigned int n_genes);
double log_pathway_prior(const vector<unsigned int> &pathway, unsigned int n_genes, unsigned int n_pathways);

// compute log likelihood for a given number of mutations given pathway size and error probabilities
double compute_log_lik_active(double n_mutations, double pathway_size, double bgp, double fbp);
double compute_log_lik_inactive(double n_mutations, double pathway_size, double bgp);

double compute_likelihood_for_sample(const vector<size_t> &r,
                                     vector<size_t> &pathway_sizes,
                                     size_t stage,
                                     double bgp, double fbp);
double compute_likelihood_for_sample(const vector<size_t> &r,
                                     const LinearProgressionState &state,
                                     size_t stage,
                                     double bgp, double fbp);

// first one is called by parameter proposal for single state
// R: number of mutations for each patient for each pathway
double compute_pathway_likelihood(vector<vector<size_t> > &R,
                                  const LinearProgressionState &state,
                                  const LinearProgressionParameters &params);
// second one is called when computing likelihood of a particle
double compute_pathway_likelihood(const gsl_matrix &obs,
                                  const vector<size_t> &sum_vec,
                                  const LinearProgressionState &state,
                                  const LinearProgressionParameters &params);
// third one is same as the second one except that it uses likelihood table cached by LPMModel.
//double compute_pathway_likelihood(gsl_matrix &obs_matrix,
//                                  vector<size_t> &sum_vec,
//                                  const LinearProgressionState &state,
//                                  const LinearProgressionParameters &params,
//                                  unordered_map<size_t, unordered_map<size_t, double> > _dict_active_probs,
//                                  unordered_map<size_t, unordered_map<size_t, double> > _dict_inactive_probs);

void fast_matrix_product(gsl_matrix &obs, size_t m, const vector<size_t> &pathway_membership, vector<size_t> &ret);

//double compute_pathway_likelihood(vector<vector<size_t> > &R, const LinearProgressionState &state, const LinearProgressionParameters &params);
//double compute_pathway_likelihood_cache_version1(gsl_matrix &obs_matrix, vector<size_t> &sum_vec, const LinearProgressionState &state, const LinearProgressionParameters &params);
//double compute_pathway_likelihood_cache_version2(gsl_matrix &obs_matrix, vector<size_t> &sum_vec, const LinearProgressionState &state, const LinearProgressionParameters &params, unordered_map<size_t, unordered_map<size_t, double>> &_dict_active_probs, unordered_map<size_t, unordered_map<size_t, double>> &_dict_inactive_probs);

//double compute_pathway_likelihood(gsl_matrix_view &Y, gsl_matrix_view &X, vector<size_t> &stages, double fbp, double bgp);
// matrix product between m-th row of obs and mutation matrix induced by pathway_membership

#endif /* lpm_likelihood_hpp */
