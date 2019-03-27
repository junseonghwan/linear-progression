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
    vector<size_t> pathway_membership; // map from gene id to pathway membership (essentially)
    vector<unordered_set<size_t> > pathways; // pathway idx to set of genes

    // MxK
    // element (i, j): number of mutations for patient in pathway j for the current state
    vector<vector<size_t> > _cache_counts;

    size_t n_genes = 0;
    size_t num_driver_pathways = 0;
    bool allocate_passenger_pathway = false;

    double log_lik = 0.0;
    
    gsl_matrix &obs;
    vector<size_t> &row_sum;
    void initialize_pathways();
    void update_cache(size_t g, size_t old_pathway, size_t new_pathway);

public:
    //LinearProgressionState(size_t n_genes, size_t num_driver_pathways, bool allocate_passenger_pathway);
    LinearProgressionState(gsl_matrix &obs, vector<size_t> &row_sums, size_t n_genes, size_t num_driver_pathways, bool allocate_passenger_pathway);
    void sample_from_prior(gsl_rng *random);
    void sample_min_valid_pathway(gsl_rng *random);
    void swap_pathways(size_t pathway1, size_t pathway2);
    void update_pathway_membership(size_t gene_idx, size_t new_pathway);
    
    void compute_counts_for_sample(gsl_matrix &obs, vector<size_t> &row_sums, size_t m, vector<size_t> &ret) const;
    size_t get_pathway_size(size_t k) const;
    size_t get_num_pathways() const;
    inline double get_log_lik() const { return log_lik; }
    inline size_t get_num_driver_pathways () const { return num_driver_pathways; }
    size_t get_pathway_membership_of(size_t gene_idx) const;
    const vector<size_t> &get_pathway_membership() const;
    bool contains_empty_driver_pathway() const;
    bool has_passenger_pathway() const;
    const vector<size_t> &get_cache_at(size_t m) const;

    inline void set_num_driver_pathways(size_t K) { this->num_driver_pathways = K; }
    inline void set_log_lik(double log_lik) { this->log_lik = log_lik; }

    LinearProgressionState *increase_num_drivers(gsl_rng *random);
    string to_string() const;

    const static size_t PASSENGER_PATHWAY_IDX = 0;
    ~LinearProgressionState();
};

#endif /* lpm_state_hpp */
