//
//  lpm_state.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#include <algorithm>
#include <iostream>
#include <iterator>
#include <gsl/gsl_randist.h>
#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "data_util.hpp"
#include "lpm_likelihood.hpp"
#include "lpm_state.hpp"

LinearProgressionState::LinearProgressionState(const gsl_matrix &obs, const vector<unsigned int> &row_sums, unsigned int n_genes, unsigned int num_driver_pathways, bool allocate_passenger_pathway) :
n_genes(n_genes), num_driver_pathways(num_driver_pathways), allocate_passenger_pathway(allocate_passenger_pathway), obs(obs), row_sum(row_sums)
{
    initialize_pathways();

    unsigned int n_patients = obs.size1;
    unsigned int n_pathways = allocate_passenger_pathway ? num_driver_pathways + 1 : num_driver_pathways;
    // initialize cache
    for (unsigned int m = 0; m < n_patients; m++) {
        _cache_counts.push_back(vector<unsigned int>(n_pathways));
        compute_counts_for_sample(obs, row_sums, m, _cache_counts[m]);
    }
}

void LinearProgressionState::initialize_pathways()
{
    // initialize pathways
    for (unsigned int k = 0; k < get_num_pathways(); k++) {
        pathways.push_back(unordered_set<unsigned int>());
    }
    // initialize pathway membership vector (to 0)
    for (unsigned int n = 0; n < n_genes; n++) {
        pathway_membership.push_back(0);
        pathways[0].insert(n);
    }
}

void LinearProgressionState::sample_min_valid_pathway(gsl_rng *random)
{
    if (!allocate_passenger_pathway) {
        // no passenger pathway, just sample from prior
        sample_pathway(random);
        return;
    }

    unsigned int *genes = new unsigned int[n_genes];
    for (unsigned int n = 0; n < n_genes; n++) {
        genes[n] = n;
    }
    // populate each of the driver pathways with exactly one gene, and allocate the remaining to passenger pathway
    gsl_ran_shuffle(random, genes, n_genes, sizeof(unsigned int));
    for (unsigned int n = 0; n < n_genes; n++) {
        if (n < num_driver_pathways) {
            update_pathway_membership(genes[n], n);
        } else {
            update_pathway_membership(genes[n], num_driver_pathways);
        }
    }
    delete [] genes;
}

void LinearProgressionState::sample_pathway(gsl_rng *random)
{
    vector<unsigned int> new_pathway_membership(n_genes);
    propose_pathway(random, get_num_pathways(), new_pathway_membership);

    for (unsigned int g = 0; g < n_genes; g++) {
        update_pathway_membership(g, new_pathway_membership[g]);
    }
}

void LinearProgressionState::update_cache(unsigned int g, unsigned int old_pathway, unsigned int new_pathway)
{
    unsigned int n_patients = obs.size1;
    for (unsigned int m = 0; m < n_patients; m++) {
        vector<unsigned int> &r = _cache_counts[m];
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

void LinearProgressionState::update_pathway_membership(unsigned int gene_idx, unsigned int new_pathway)
{
    unsigned int old_pathway = pathway_membership[gene_idx];
    pathway_membership[gene_idx] = new_pathway;
    pathways[old_pathway].erase(gene_idx);
    pathways[new_pathway].insert(gene_idx);

    // update the cache
    if (old_pathway != new_pathway) {
        update_cache(gene_idx, old_pathway, new_pathway);
    }
}

void LinearProgressionState::swap_pathways(unsigned int i, unsigned int j)
{
    if (i > get_num_pathways() || j > get_num_pathways()) {
        cerr << i << " or " << j << " > " << get_num_pathways() << endl;
        exit(-1);
    }
    unordered_set<unsigned int> temp_i(pathways[i]);
    unordered_set<unsigned int> temp_j(pathways[j]);
    for (auto it = temp_i.begin(); it != temp_i.end(); ++it) {
        update_pathway_membership(*it, j);
    }
    for (auto it = temp_j.begin(); it != temp_j.end(); ++it) {
        update_pathway_membership(*it, i);
    }
}

unsigned int LinearProgressionState::get_pathway_size(unsigned int k) const
{
    if (k >= get_num_pathways()) {
        cerr << "Error: Index out of bounds. k: " << k << " num_pathways: " << get_num_pathways() << endl;
        exit(-1);
    }
    return pathways[k].size();
}

unsigned int LinearProgressionState::get_num_pathways() const
{
    if (allocate_passenger_pathway) {
        return num_driver_pathways + 1;
    }
    return num_driver_pathways;
}

unsigned int LinearProgressionState::get_n_patients() const
{
    return obs.size1;
}

bool LinearProgressionState::has_passenger_pathway() const
{
    return allocate_passenger_pathway;
}

bool LinearProgressionState::contains_empty_driver_pathway() const
{
    for (unsigned int k = 0; k < num_driver_pathways; k++) {
        if (pathways[k].size() == 0) {
            return true;
        }
    }
    return false;
}

void LinearProgressionState::compute_counts_for_sample(const gsl_matrix &obs_matrix, const vector<unsigned int> &row_sums, unsigned int m, vector<unsigned int> &r) const
{
    if (r.size() != get_num_pathways()) {
        cerr << "Error: return vector size does not match the number of pathways." << endl;
        exit(-1);
    }
    unsigned int sum = 0;
    for (unsigned int k = 0; k < num_driver_pathways; k++) {
        r[k] = 0; // reset the counts
        for (unsigned int col_idx : pathways[k]) {
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
    for (unsigned int n = 0; n < n_genes; n++) {
        unsigned int pathway = this->get_pathway_membership_of(n);
        new_state->update_pathway_membership(n, pathway);
    }

    // if the passenger pathway of this object is empty, then we need to sample a gene
    // move it to the new pathway as long as moving that gene does not lead to creating an empty pathway
    if (new_state->get_pathway_size(this->num_driver_pathways) == 0) {
        vector<double> probs(this->num_driver_pathways, 0.0); // from one of new_K - 1 pathways, move it to the new driver pathway
        double sum = 0;
        for (unsigned int k = 0; k < this->num_driver_pathways; k++) {
            if (new_state->get_pathway_size(k) >= 2) {
                probs[k] = 1.0;
                sum++;
            }
        }
        for (unsigned int k = 0; k < this->num_driver_pathways; k++) {
            probs[k] /= sum;
        }
        unsigned int k = multinomial(random, probs);
        double u = gsl_ran_flat(random, 0.0, 1.0);
        double incr = 1.0/pathways[k].size();
        sum = 0.0;
        for (unsigned int g : new_state->pathways[k]) {
            if (u < sum + incr) {
                new_state->update_pathway_membership(g, this->num_driver_pathways);
                break;
            }
            sum += incr;
        }
    }
    return new_state;
}

unsigned int LinearProgressionState::get_pathway_membership_of(unsigned int gene_idx) const
{
    return pathway_membership[gene_idx];
}

const vector<unsigned int> &LinearProgressionState::get_pathway_membership() const
{
    return pathway_membership;
}

string LinearProgressionState::to_string() const
{
    string str = "";
    for (unsigned int i = 0; i < n_genes; i++) {
        str += std::to_string(pathway_membership[i]);
        if (i < (n_genes - 1))
             str += ", ";
    }
    return str;
}

const vector<unsigned int> &LinearProgressionState::get_cache_at(unsigned int m) const
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

