//
//  lpm_likelihood.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-08.
//

#include <cmath>

#include <gsl/gsl_sf.h>

#include <spf/numerical_utils.hpp>

#include "lpm_likelihood.hpp"

double log_pathway_uniform_prior(unsigned int n_pathways, unsigned int n_genes)
{
    vector<double> f(n_pathways);
    f[0] = 1.0;
    for (unsigned int k = 1; k < n_pathways; k++) {
        f[k] = pow(k+1, n_genes);
        for (unsigned int j = 0; j < k; j++) {
            f[k] -= gsl_sf_choose(k+1, j+1) * f[j];
        }
    }
    return -log(f[n_pathways-1]);
}

double log_min_valid_pathway_proposal(const LinearProgressionState &pathway, unsigned int n_genes)
{
    if (pathway.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }
    double log_proposal = 0.0;
    for (unsigned int k = 0; k < pathway.get_num_driver_pathways(); k++) {
        log_proposal += log(n_genes - k);
    }
    return -log_proposal;
}


double log_pathway_proposal(const LinearProgressionState &pathway, unsigned int n_genes)
{
    if (pathway.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }
    unsigned int n_pathways = pathway.get_num_pathways();
    double log_prior = 0.0;
    for (unsigned int k = 0; k < pathway.get_num_pathways(); k++) {
        log_prior += gsl_sf_lnfact(pathway.get_pathway_size(k));
    }
    log_prior -= gsl_sf_lnchoose(n_genes - 1, n_pathways - 1);
    log_prior -= gsl_sf_lnfact(n_genes);
    return log_prior;
}

double log_pathway_proposal(const vector<unsigned int> &pathway, unsigned int n_pathways)
{
    unsigned int n_genes = pathway.size();
    unsigned int *pathway_sizes = new unsigned int[n_pathways];
    for (unsigned int k = 0; k < n_pathways; k++)
        pathway_sizes[k] = 0;
    for (unsigned int g = 0; g < n_genes; g++) {
        pathway_sizes[pathway[g]] += 1;
    }
    
    double log_prior = 0.0;
    for (unsigned int k = 0; k < n_pathways; k++) {
        log_prior += gsl_sf_lnfact(pathway_sizes[k]);
    }
    log_prior -= gsl_sf_lnchoose(n_genes - 1, n_pathways - 1);
    log_prior -= gsl_sf_lnfact(n_genes);
    
    delete [] pathway_sizes;

    return log_prior;
}


double compute_log_lik_active(unsigned int n_mutations, unsigned int pathway_size, double bgp, double fbp)
{
    if (pathway_size <= 0) {
        // this is not in the support set
        return DOUBLE_NEG_INF;
    }
    double log_lik = 0.0;
    if (n_mutations == 0) {
        // one entry got flipped back from 1 -> 0 and the rest are not impacted by back ground prob
        log_lik += log(fbp);
        log_lik += (pathway_size - 1) * log(1-bgp);
        log_lik += log(pathway_size);
    } else if (n_mutations == 1) {
        // no error case: one entry was unaffected by fbp (1->1) and the rest are not affected by background prob
        log_lik += log(1 - fbp);
        log_lik += (pathway_size - 1) * log(1 - bgp);
        if (pathway_size > n_mutations) {
            // we need to account for flip back (1->0) of a gene
            // followed by background mutation of the one that gave rise to n_mutations = 1
            // remaining genes are unaffected
            // e.g., obs may be 1,0,0, which may have been obtained by
            // fb and bg mutation from 0,1,0 or 0,0,1
            double log_lik2 = log(fbp) + log(bgp);
            log_lik2 += (pathway_size - 2) * log(1 - bgp);
            log_lik2 += log(pathway_size - 1);
            log_lik = log_add(log_lik, log_lik2);
        }
    } else {
        // one entry was unaffected by fbp (1->1)
        // (n_mutations - 1) are background mutations
        // the rest are not affected by background prob
        // e.g., 1,1,1,0 can be obtained from 1,0,0,0 with bgm on second and third entry
        // but we do not know which one is the unaffected one, so the number of possibilities are
        // (1,0,0,0), (0,1,0,0) or (0,0,1,0) with 2 bg mutations, no flip back, and no bg on the last entry.
        log_lik += log(n_mutations); // multiply by the number of possibilities
        log_lik += log(1 - fbp);
        log_lik += (n_mutations - 1) * log(bgp);
        log_lik += (pathway_size - n_mutations) * log(1 - bgp);
        if (pathway_size > n_mutations) {
            // account for flip back (1->0) followed by n_mutations bgm
            // e.g., 0,0,0,1->0,0,0,0->1,1,1,0
            // note: we don't allow fb followed by bgm
            double log_lik2 = log(pathway_size - n_mutations);
            log_lik2 += log(fbp);
            log_lik2 += n_mutations * log(bgp);
            log_lik2 += (pathway_size - n_mutations - 1) * log(1 - bgp);
            log_lik = log_add(log_lik, log_lik2);
        }
    }
    log_lik -= log(pathway_size);
    return log_lik;
}

double compute_log_lik_inactive(unsigned int n_mutations, unsigned int pathway_size, double bgp)
{
    // all mutations are background mutations
    double log_lik = n_mutations * log(bgp);
    log_lik += (pathway_size - n_mutations) * log(1 - bgp);
    return log_lik;
}

double compute_likelihood_for_sample(const vector<unsigned short> &r, const LinearProgressionState &state, unsigned int stage, double bgp, double fbp)
{
    double log_lik = 0.0;
    unsigned int K = state.get_num_pathways();
    for (int k = 0; k < K; k++) {
        double val = r[k];
        if (k <= stage) {
            log_lik += compute_log_lik_active(val, state.get_pathway_size(k), bgp, fbp);
        } else {
            log_lik += compute_log_lik_inactive(val, state.get_pathway_size(k), bgp);
        }
    }
    return log_lik;
}

void reset_r_vector(vector<unsigned int> &vec)
{
    for (unsigned int i = 0; i < vec.size(); i++) {
        vec[i] = 0;
    }
}

bool contains_empty_pathway(vector<unsigned int> &pathway_sizes)
{
    for (unsigned int i = 0; i < pathway_sizes.size(); i++) {
        if (pathway_sizes[i] == 0) {
            return true;
        }
    }
    return false;
}

double compute_pathway_likelihood(const LinearProgressionState &state,
                                  const LinearProgressionParameters &params)
{
    double log_lik = 0.0;

    unsigned int n_obs = state.get_n_patients();
    unsigned int n_driver_pathways = state.get_num_driver_pathways();
    double log_n_driver_pathways = log(n_driver_pathways);
    if (state.contains_empty_driver_pathway()) {
        return DOUBLE_NEG_INF;
    }

    //marginalize out stages for each patient
    for (unsigned int m = 0; m < n_obs; m++) {
        const vector<unsigned short> &ret = state.get_cache_at(m);
        double log_lik_stage = compute_likelihood_for_sample(ret, state, 0, params.get_bgp(), params.get_fbp());
        double log_lik_m = log_lik_stage - log_n_driver_pathways; // -log_n_driver_pathways accounts for prior on stage assignment: (1/n_pathways)
        // marginalize over patient stages taking on {1, ..., driver pathways}
        for (unsigned int stage = 1; stage < state.get_num_driver_pathways(); stage++) {
            unsigned int pathway_size = state.get_pathway_size(stage);
            unsigned int r = ret[stage];
            log_lik_stage += compute_log_lik_active(r, pathway_size, params.get_bgp(), params.get_fbp());
            log_lik_stage -= compute_log_lik_inactive(r, pathway_size, params.get_bgp());
            log_lik_m = log_add(log_lik_m, log_lik_stage - log_n_driver_pathways);
        }
        log_lik += log_lik_m;
    }
    return log_lik;
}

void fast_matrix_product(gsl_matrix &obs, unsigned int m, const vector<unsigned int> &pathway_membership, vector<unsigned int> &r)
{
    unsigned int n_genes = pathway_membership.size();
    for (unsigned int n = 0; n < n_genes; n++) {
        double val = gsl_matrix_get(&obs, m, n);
        if (val == 1.0) {
            unsigned int pathway_idx = pathway_membership[n];
            if (pathway_idx >= r.size()) {
                cerr << "Error in fast_matrix_product: invalid size of ret vector is allocated." << endl;
                exit(-1);
            }
            r[pathway_idx] += 1;
        }
    }
}


