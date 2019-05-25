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
    vector<unsigned int> pathway_membership; // map from gene id to pathway membership (essentially)
    vector<unordered_set<unsigned int> > pathways; // pathway idx to set of genes

    // MxK
    // element (i, j): number of mutations for patient i in pathway j for the current state
    // this cache gets updated when gene membership changes
    vector<vector<unsigned int> > _cache_counts;

    unsigned int n_genes = 0;
    unsigned int num_driver_pathways = 0;
    bool allocate_passenger_pathway = false;

    double log_lik = 0.0;

    void initialize_pathways();

    const gsl_matrix &obs;
    //const vector<unsigned int> &row_sum;
    void update_cache(unsigned int g, unsigned int old_pathway, unsigned int new_pathway);

public:
    LinearProgressionState(const gsl_matrix &obs, const vector<unsigned int> &row_sums, unsigned int n_genes, unsigned int num_driver_pathways, bool allocate_passenger_pathway);
    void sample_pathway(gsl_rng *random);
    void sample_min_valid_pathway(gsl_rng *random);
    void swap_pathways(unsigned int pathway1, unsigned int pathway2);
    void update_pathway_membership(unsigned int gene_idx, unsigned int new_pathway);
    
    void compute_counts_for_sample(const gsl_matrix &obs, const vector<unsigned int> &row_sums, unsigned int m, vector<unsigned int> &ret) const;
    unsigned int get_pathway_size(unsigned int k) const;
    unsigned int get_num_pathways() const;
    unsigned int get_n_patients() const;
    inline double get_log_lik() const { return log_lik; }
    inline unsigned int get_num_driver_pathways () const { return num_driver_pathways; }
    unsigned int get_pathway_membership_of(unsigned int gene_idx) const;
    const vector<unsigned int> &get_pathway_membership() const;
    bool contains_empty_driver_pathway() const;
    bool has_passenger_pathway() const;
    const vector<unsigned int> &get_cache_at(unsigned int m) const;

    inline void set_num_driver_pathways(unsigned int K) { this->num_driver_pathways = K; }
    inline void set_log_lik(double log_lik) { this->log_lik = log_lik; }

    //LinearProgressionState *increase_num_drivers(gsl_rng *random);
    string to_string() const;

    inline bool operator==(const LinearProgressionState& other) const
    {
        return (this->to_string() == other.to_string());
    }
    
    ~LinearProgressionState();
};

namespace std {
    template<>
    struct hash<LinearProgressionState>
    {
        std::size_t operator()(const LinearProgressionState& key) const
        {
            size_t hash_key = hash<string>()(key.to_string());
            return hash_key;
        }
    };
}


#endif /* lpm_state_hpp */
