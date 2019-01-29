//
//  lpm_state.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#include <gsl/gsl_randist.h>

#include "lpm_state.hpp"

LinearProgressionState::LinearProgressionState(size_t n_genes, size_t num_pathways, bool allocate_passenger_pathway) :
n_genes(n_genes), num_pathways(num_pathways), allocate_passenger_pathway(allocate_passenger_pathway)
{
    initialize_pathway_membership();
}

void LinearProgressionState::initialize_pathway_membership()
{
    for (size_t k = 0; k < num_pathways; k++) {
        pathway_sizes.push_back(0);
    }
    for (size_t n = 0; n < n_genes; n++) {
        pathway_membership.push_back(0);
    }
    pathway_sizes[0] = n_genes;
}

void LinearProgressionState::sample_from_prior(gsl_rng *random)
{
    size_t *genes = new size_t[n_genes];
    for (size_t n = 0; n < n_genes; n++) {
        genes[n] = n;
    }

    // there must be at least one gene in each pathway (empty pathway has zero measure)
    // 1. shuffle the genes
    // 2. first K = num_pathway genes are allocated to each of K pathways
    // 3. the remaining genes sample a pathway at random
    gsl_ran_shuffle(random, genes, n_genes, sizeof(size_t));
    for (size_t n = 0; n < n_genes; n++) {
        size_t g = genes[n];
        if (n < num_pathways) {
            update_pathway_membership(g, n);
        } else {
            size_t k = gsl_rng_uniform_int(random, num_pathways);
            update_pathway_membership(g, k);
        }
    }
}

void LinearProgressionState::update_pathway_membership(size_t gene_idx, size_t pathway)
{
    size_t old_pathway = pathway_membership[gene_idx];
    pathway_membership[gene_idx] = pathway;
    pathway_sizes[old_pathway] -= 1;
    pathway_sizes[pathway] += 1;
}

vector<size_t> LinearProgressionState::get_pathway()
{
    return pathway_membership;
}

vector<size_t> LinearProgressionState::get_pathway_sizes()
{
    return pathway_sizes;
}

string LinearProgressionState::to_string()
{
    string str = "";
    for (size_t i = 0; i < n_genes; i++) {
        str += std::to_string(pathway_membership[i]);
        if (i < (n_genes - 1))
             str += ", ";
    }
    return str;
}
