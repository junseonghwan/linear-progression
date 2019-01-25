//
//  lpm_state.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#include "lpm_state.hpp"

LinearProgressionState::LinearProgressionState(size_t n_genes, size_t num_pathways, bool allocate_passenger_pathway) :
n_genes(n_genes), num_pathways(num_pathways), allocate_passenger_pathway(allocate_passenger_pathway)
{
}

LinearProgressionState *LinearProgressionState::copy_state()
{
    LinearProgressionState *copy_state = new LinearProgressionState(this->n_genes, this->num_pathways, this->allocate_passenger_pathway);
    copy_state->pathway_membership = this->pathway_membership;
    copy_state->log_lik = this->get_log_lik();
    return copy_state;
}

//LinearProgressionState &LinearProgressionState::operator= ( LinearProgressionState & other)
//{
//    this->num_pathways = other.num_pathways;
//    this->n_genes = other.n_genes;
//    this->pathway_membership = other.pathway_membership;
//    this->pathway_membership[0] = 2;
//    this->log_lik = other.log_lik;
//    return *this;
//}

void LinearProgressionState::initialize(gsl_rng *random)
{
    for (size_t n = 0; n < n_genes; n++) {
        // sample pathway membership
        size_t k = gsl_rng_uniform_int(random, num_pathways);
        pathway_membership.push_back(k);
    }
}
void LinearProgressionState::update_pathway_membership(size_t gene_idx, size_t pathway)
{
    pathway_membership[gene_idx] =  pathway;
}

vector<size_t> &LinearProgressionState::get_pathway()
{
    return pathway_membership;
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
