//
//  lpm_state.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#include <iostream>
#include <gsl/gsl_randist.h>

#include "lpm_likelihood.hpp"
#include "lpm_state.hpp"

LinearProgressionState::LinearProgressionState(size_t n_genes, size_t num_driver_pathways, bool allocate_passenger_pathway) :
n_genes(n_genes), num_driver_pathways(num_driver_pathways), allocate_passenger_pathway(allocate_passenger_pathway)
{
    unordered_set<size_t> *pathway = 0;
    for (size_t k = 0; k < get_num_pathways(); k++) {
        pathway = new unordered_set<size_t>();
        pathways.push_back(*pathway);
    }
    for (size_t n = 0; n < n_genes; n++) {
        pathway_membership.push_back(0);
        pathways[0].insert(n);
    }
}

void LinearProgressionState::sample_min_valid_pathway(gsl_rng *random)
{
    if (!allocate_passenger_pathway) {
        // no passenger pathway, just sample from prior
        sample_from_prior(random);
        return;
    }
    size_t *genes = new size_t[n_genes];
    for (size_t n = 0; n < n_genes; n++) {
        genes[n] = n;
    }

    // populate each of the driver pathways with exactly one gene, and allocate the remaining to passenger pathway
    gsl_ran_shuffle(random, genes, n_genes, sizeof(size_t));
    for (size_t n = 0; n < n_genes; n++) {
        if (n < num_driver_pathways) {
            update_pathway_membership(genes[n], n);
        } else {
            update_pathway_membership(genes[n], num_driver_pathways);
        }
    }
}

void LinearProgressionState::sample_from_prior(gsl_rng *random)
{
    size_t *genes = new size_t[n_genes];
    for (size_t n = 0; n < n_genes; n++) {
        genes[n] = n;
    }

    // ensure there to be at least one gene in each pathway (empty pathway has zero measure)
    // 1. shuffle the genes
    // 2. place the first num_driver_pathway genes to first num_driver_pathway pathways
    // 3. place remaining genes randomly (including passenger pathway)
    size_t num_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
    gsl_ran_shuffle(random, genes, n_genes, sizeof(size_t));
    for (size_t n = 0; n < n_genes; n++) {
        size_t g = genes[n];
        if (n < num_driver_pathways) {
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
    pathways[old_pathway].erase(gene_idx);
    pathways[pathway].insert(gene_idx);
}

//vector<size_t> LinearProgressionState::get_pathway()
//{
//    return pathway_membership;
//}
//
//vector<size_t> LinearProgressionState::get_pathway_sizes()
//{
//    return pathway_sizes;
//}

void LinearProgressionState::swap_pathways(size_t i, size_t j)
{
    if (i > get_num_pathways() || j > get_num_pathways()) {
        cerr << i << " or " << j << " > " << get_num_pathways() << endl;
        exit(-1);
    }
    for (size_t g1 : pathways[i]) {
        pathway_membership[g1] = j;
    }
    for (size_t g2 : pathways[j]) {
        pathway_membership[g2] = i;
    }
    unordered_set<size_t> temp = pathways[j];
    pathways[j] = pathways[i];
    pathways[i] = temp;
}

size_t LinearProgressionState::get_pathway_size(size_t k)
{
    if (k >= get_num_pathways()) {
        cerr << "Error: Index out of bounds. k: " << k << " num_pathways: " << get_num_pathways() << endl;
        exit(-1);
    }
    return pathways[k].size();
}

size_t LinearProgressionState::get_num_pathways()
{
    if (allocate_passenger_pathway) {
        return num_driver_pathways + 1;
    }
    return num_driver_pathways;
}

bool LinearProgressionState::has_passenger_pathway()
{
    return allocate_passenger_pathway;
}

bool LinearProgressionState::contains_empty_driver_pathway()
{
    for (size_t k = 0; k < num_driver_pathways; k++) {
        if (pathways[k].size() == 0) {
            return true;
        }
    }
    return false;
}

void LinearProgressionState::compute_counts_for_sample(gsl_matrix_view &obs, vector<size_t> &row_sums, size_t m, vector<size_t> &r)
{
    if (r.size() != get_num_pathways()) {
        cerr << "Error: return vector size does not match the number of pathways." << endl;
        exit(-1);
    }
    size_t sum = 0;
    for (size_t k = 0; k < num_driver_pathways; k++) {
        r[k] = 0; // reset the counts
        for (size_t col_idx : pathways[k]) {
            double val = gsl_matrix_get(&obs.matrix, m, col_idx);
            if (val == 1.0) {
                r[k] += 1;
            }
        }
        sum += r[k];
    }
    if (allocate_passenger_pathway) {
        r[num_driver_pathways] = row_sums[m] - sum;
    } else {
        // sanity check
        if (row_sums[m] != sum) {
            cerr << "Error in counting number of mutations for sample " << m << endl;
            exit(-1);
        }
    }
}

LinearProgressionState *LinearProgressionState::LinearProgressionState::increase_num_drivers()
{
    LinearProgressionState *new_state = new LinearProgressionState(this->n_genes, this->num_driver_pathways + 1, this->allocate_passenger_pathway);
    for (size_t n = 0; n < n_genes; n++) {
        size_t pathway = this->get_pathway_membership(n);
        new_state->update_pathway_membership(n, pathway);
    }    
    return new_state;
}

size_t LinearProgressionState::get_pathway_membership(size_t gene_idx)
{
    return pathway_membership[gene_idx];
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
