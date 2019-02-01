//
//  lpm_state.hpp
//
//  Maintains a map from gene id to pathway membership backed by vector.
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#ifndef lpm_state_hpp
#define lpm_state_hpp

#include <stdio.h>
#include <string>
#include <unordered_set>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

using namespace std;

class LinearProgressionState
{
    //gsl_matrix *pathways;
    vector<size_t> pathway_membership; // map from gene id to pathway membership (essentially)
    vector<unordered_set<size_t>> pathways; // pathway to set of genes
    //vector<size_t> pathway_sizes;

    size_t n_genes = 0;
    size_t num_driver_pathways = 0;
    bool allocate_passenger_pathway = false;

    double log_lik = 0.0;

public:
    LinearProgressionState(size_t n_genes, size_t num_driver_pathways, bool allocate_passenger_pathway = false);
    void sample_from_prior(gsl_rng *random);
    void sample_min_valid_pathway(gsl_rng *random);
    void swap_pathways(size_t pathway1, size_t pathway2);
    void update_pathway_membership(size_t gene_idx, size_t pathway);
    void compute_counts_for_sample(gsl_matrix_view &obs, vector<size_t> &row_sums, size_t m, vector<size_t> &ret);
    //vector<size_t> get_pathway();
    //vector<size_t> get_pathway_sizes();
    size_t get_pathway_size(size_t k);

    size_t get_num_pathways();
    inline size_t get_num_driver_pathways() { return num_driver_pathways; }
    inline void set_num_driver_pathways(size_t K) { this->num_driver_pathways = K; }
    inline void set_log_lik(double log_lik) { this->log_lik = log_lik; }
    inline double get_log_lik() { return log_lik; }
    bool has_passenger_pathway();
    bool contains_empty_driver_pathway();
    LinearProgressionState *increase_num_drivers();
    size_t get_pathway_membership(size_t gene_idx);
    string to_string();
    
    const static size_t PASSENGER_PATHWAY_IDX = 0;
};

#endif /* lpm_state_hpp */
