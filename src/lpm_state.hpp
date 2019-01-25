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
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

using namespace std;

class LinearProgressionState
{
    //gsl_matrix *pathways;
    vector<size_t> pathway_membership; // map from gene id to pathway membership (essentially)

    size_t n_genes = 0;
    size_t num_pathways = 0;
    bool allocate_passenger_pathway = false;
    
    double log_lik = 0.0;
    
public:
    //inline LinearProgressionState() { }
    LinearProgressionState(size_t n_genes, size_t num_pathways, bool allocate_passenger_pathway = false);
    //LinearProgressionState &operator= ( LinearProgressionState & );
    LinearProgressionState *copy_state();
    void initialize(gsl_rng *random);
    void update_pathway_membership(size_t gene_idx, size_t pathway);
    vector<size_t> &get_pathway();

    inline size_t get_num_pathways() { return num_pathways; }
    inline void set_num_pathways(size_t K) { this->num_pathways = K; }
    inline void set_log_lik(double log_lik) { this->log_lik = log_lik; }
    inline double get_log_lik() { return log_lik; }
    string to_string();
    
    const size_t PASSENGER_PATHWAY_IDX = 0;
};

#endif /* lpm_state_hpp */
