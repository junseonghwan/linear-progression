//
//  lpm_state.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#include <iostream>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_state.hpp"

//LinearProgressionState::LinearProgressionState(size_t n_genes, size_t num_driver_pathways, bool allocate_passenger_pathway) :
//n_genes(n_genes), num_driver_pathways(num_driver_pathways), allocate_passenger_pathway(allocate_passenger_pathway)
//{
//    initialize_pathways();
//}

void LinearProgressionState::initialize_pathways()
{
    // initialize pathways
    for (size_t k = 0; k < get_num_pathways(); k++) {
        pathways.push_back(unordered_set<size_t>());
    }
    // initialize pathway membership vector (to 0)
    for (size_t n = 0; n < n_genes; n++) {
        pathway_membership.push_back(0);
        pathways[0].insert(n);
    }
}

LinearProgressionState::LinearProgressionState(gsl_matrix &obs, vector<size_t> &row_sums, size_t n_genes, size_t num_driver_pathways, bool allocate_passenger_pathway) :
n_genes(n_genes), num_driver_pathways(num_driver_pathways), allocate_passenger_pathway(allocate_passenger_pathway), obs(obs), row_sum(row_sums)
{
    initialize_pathways();

    size_t n_patients = obs.size1;
    size_t n_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
    // initialize cache
    for (size_t m = 0; m < n_patients; m++) {
        _cache_counts.push_back(vector<size_t>(n_pathways));
        compute_counts_for_sample(obs, row_sums, m, _cache_counts[m]);
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
    delete [] genes;
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
    
    delete [] genes;
}

void LinearProgressionState::update_cache(size_t g, size_t old_pathway, size_t new_pathway)
{
    size_t n_patients = obs.size1;
    for (size_t m = 0; m < n_patients; m++) {
        vector<size_t> &r = _cache_counts[m];
        double val = gsl_matrix_get(&obs, m, g);
        if (val == 1.0) {
            if (r[old_pathway] <= 0) {
                cerr << "Error: cache computation has a bug." << endl;
                exit(-1);
            }
            r[old_pathway] -= 1;
            r[new_pathway] += 1;
        }
    }
}

void LinearProgressionState::update_pathway_membership(size_t gene_idx, size_t new_pathway)
{
    size_t old_pathway = pathway_membership[gene_idx];
    pathway_membership[gene_idx] = new_pathway;
    pathways[old_pathway].erase(gene_idx);
    pathways[new_pathway].insert(gene_idx);

    // update the cache
    if (old_pathway != new_pathway) {
        update_cache(gene_idx, old_pathway, new_pathway);
    }
}

void LinearProgressionState::swap_pathways(size_t i, size_t j)
{
    if (i > get_num_pathways() || j > get_num_pathways()) {
        cerr << i << " or " << j << " > " << get_num_pathways() << endl;
        exit(-1);
    }
    unordered_set<size_t> temp_i(pathways[i]);
    unordered_set<size_t> temp_j(pathways[j]);
    for (auto it = temp_i.begin(); it != temp_i.end(); ++it) {
        update_pathway_membership(*it, j);
    }
    for (auto it = temp_j.begin(); it != temp_j.end(); ++it) {
        update_pathway_membership(*it, i);
    }
}

size_t LinearProgressionState::get_pathway_size(size_t k) const
{
    if (k >= get_num_pathways()) {
        cerr << "Error: Index out of bounds. k: " << k << " num_pathways: " << get_num_pathways() << endl;
        exit(-1);
    }
    return pathways[k].size();
}

size_t LinearProgressionState::get_num_pathways() const
{
    if (allocate_passenger_pathway) {
        return num_driver_pathways + 1;
    }
    return num_driver_pathways;
}

bool LinearProgressionState::has_passenger_pathway() const
{
    return allocate_passenger_pathway;
}

bool LinearProgressionState::contains_empty_driver_pathway() const
{
    for (size_t k = 0; k < num_driver_pathways; k++) {
        if (pathways[k].size() == 0) {
            return true;
        }
    }
    return false;
}

void LinearProgressionState::compute_counts_for_sample(gsl_matrix &obs_matrix, vector<size_t> &row_sums, size_t m, vector<size_t> &r) const
{
    if (r.size() != get_num_pathways()) {
        cerr << "Error: return vector size does not match the number of pathways." << endl;
        exit(-1);
    }
    size_t sum = 0;
    for (size_t k = 0; k < num_driver_pathways; k++) {
        r[k] = 0; // reset the counts
        for (size_t col_idx : pathways[k]) {
            double val = gsl_matrix_get(&obs_matrix, m, col_idx);
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

LinearProgressionState *LinearProgressionState::LinearProgressionState::increase_num_drivers(gsl_rng *random)
{
    LinearProgressionState *new_state = new LinearProgressionState(obs, row_sum, this->n_genes, this->num_driver_pathways + 1, this->allocate_passenger_pathway);
    for (size_t n = 0; n < n_genes; n++) {
        size_t pathway = this->get_pathway_membership_of(n);
        new_state->update_pathway_membership(n, pathway);
    }

    // if the passenger pathway of this object is empty, then we need to sample a gene
    // move it to the new pathway as long as moving that gene does not lead to creating an empty pathway
    if (new_state->get_pathway_size(this->num_driver_pathways) == 0) {
        vector<double> probs(this->num_driver_pathways, 0.0); // from one of new_K - 1 pathways, move it to the new driver pathway
        double sum = 0;
        for (size_t k = 0; k < this->num_driver_pathways; k++) {
            if (new_state->get_pathway_size(k) >= 2) {
                probs[k] = 1.0;
                sum++;
            }
        }
        for (size_t k = 0; k < this->num_driver_pathways; k++) {
            probs[k] /= sum;
        }
        size_t k = multinomial(random, probs);
        double u = gsl_ran_flat(random, 0.0, 1.0);
        double incr = 1.0/pathways[k].size();
        sum = 0.0;
        for (size_t g : new_state->pathways[k]) {
            if (u < sum + incr) {
                new_state->update_pathway_membership(g, this->num_driver_pathways);
                break;
            }
            sum += incr;
        }
    }
    return new_state;
}

size_t LinearProgressionState::get_pathway_membership_of(size_t gene_idx) const
{
    return pathway_membership[gene_idx];
}

const vector<size_t> &LinearProgressionState::get_pathway_membership() const
{
    return pathway_membership;
}

string LinearProgressionState::to_string() const
{
    string str = "";
    for (size_t i = 0; i < n_genes; i++) {
        str += std::to_string(pathway_membership[i]);
        if (i < (n_genes - 1))
             str += ", ";
    }
    return str;
}

const vector<size_t> &LinearProgressionState::get_cache_at(size_t m) const
{
    if (m >= _cache_counts.size()) {
        cerr << "Error: cache size is " << _cache_counts.size() << ", but accessing element at " << m << endl;
        exit(-1);
    }
    return _cache_counts[m];
}

LinearProgressionState::~LinearProgressionState()
{
}

