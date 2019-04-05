//
//  tests.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-28.
//

#include "tests.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/filesystem.hpp>

#include <spf/numerical_utils.hpp>
#include <spf/sampling_utils.hpp>

#include "data_util.hpp"
#include "lpm.hpp"
#include "lpm_likelihood.hpp"

using namespace std;

// 1. test the likelihood calculation
void test_likelihood(string data_path)
{
    string rep_path;
    string gold_standard_file_name;
    string data_file_name;
    string pathway_file_name;
    string param_file_name;
    
    unsigned int n_patients;
    unsigned int n_genes;
    unsigned int n_pathways;
    
    double gold_log_lik;
    double log_lik;
    double fbp, bgp;
    
    size_t idx;

    // free these objects within the loop
    gsl_matrix *data_matrix = 0;
    gsl_matrix *gold = 0;
    unsigned int *true_pathway = 0;

    for (size_t rep = 0; rep < 99; rep++) {
        rep_path = data_path + "/rep" + to_string(rep);
        gold_standard_file_name = rep_path + "/gold_standard.csv";
        data_file_name = rep_path + "/matrix.csv";
        pathway_file_name = rep_path + "/generative_mem_mat.csv";
        param_file_name = rep_path + "/parameters.csv";

        // read the data matrix
        data_matrix = read_data(data_file_name);
        n_genes = data_matrix->size2;

        // read the gold standard likelihood file
        gold = read_csv(gold_standard_file_name, false);

        // read the true pathway
        true_pathway = new unsigned int[n_genes];
        n_pathways = read_ground_truth_pathway_from_matrix(pathway_file_name, true_pathway);
        
        // read the parameters
        read_error_params(param_file_name, fbp, bgp);

        for (idx = 0; idx < gold->size1; idx++) {
            n_patients = gsl_matrix_get(gold, idx, 0);
            gold_log_lik = gsl_matrix_get(gold, idx, 1);

            gsl_matrix_view sub_matrix = gsl_matrix_submatrix(data_matrix, 0, 0, n_patients, n_genes);
            log_lik = compute_likelihood_from_matrix(&sub_matrix.matrix, true_pathway, n_pathways, n_genes, false, fbp, bgp);
            assert(abs(log_lik - gold_log_lik) < LIKELIHOOD_ERR_TOL);
        }
        
        // free the allocated memory
        delete [] true_pathway;
        delete data_matrix;
        delete gold;
    }
    
    cout << "Likelihood test passed!" << endl;
}

// run SMC to estimate the log marginal likelihood
void test_marginal_likelihood_estimate_smc(long seed, gsl_matrix *data_matrix, unsigned int n_pathways, unsigned int n_particles, unsigned int n_smc_iter, double fbp, double bgp, unordered_map<string, double> true_probs, double exact_log_marginal_lik)
{
    double logZ;
    double swap_prob = 0.2;
    unsigned int n_mh_iter = 5;
    unsigned int n_genes = data_matrix->size2;

    unsigned int *states = new unsigned int[n_particles * n_genes];
    double *log_weights = new double[n_particles];
    // run importance sampling with proposal from the prior
    logZ = run_smc_from_matrix(seed, data_matrix, n_pathways, n_particles, n_smc_iter, n_mh_iter, false, swap_prob, fbp, bgp, states, log_weights);
    cout << "Estimated logZ: " << logZ << endl;
    assert(abs(exact_log_marginal_lik - logZ)/abs(exact_log_marginal_lik) < SMC_MARGINAL_ERR_TOL);
    
    // perform tests using states and log_weights
    string pathway_str;
    size_t idx = 0;
    unordered_map<string, double> counts;
    normalize_destructively(log_weights, n_particles);
    for (size_t i = 0; i < n_particles; i++) {
        stringstream ss;
        for (size_t j = 0; j < n_genes; j++) {
            idx = i * n_genes + j;
            ss << states[idx];
            if (j < n_genes - 1) {
                ss << ", ";
            }
        }
        pathway_str = ss.str();
        if (!counts.count(pathway_str)) {
            counts[pathway_str] = 0;
        }
        counts[pathway_str] += log_weights[i];
        cout << pathway_str << ", " << counts[pathway_str] << endl;
    }

    double prob_est;
    double prob_truth;
    for (auto it = true_probs.begin(); it != true_probs.end(); ++it) {
        pathway_str = it->first;
        prob_est = counts[pathway_str];
        prob_truth = true_probs[pathway_str];
        cout << pathway_str << ", " << prob_est << ", " << prob_truth << endl;
        assert(abs(prob_est - prob_truth) < SMC_MARGINAL_ERR_TOL);
    }
    
    cout << "SMC test passed!" << endl;
}

void test_prior()
{
    gsl_rng *random = generate_random_object(123);

    unsigned int n_pathways = 2;
    unsigned int n_genes = 3;
    size_t n_mc_samples = 100000;
    
    vector<unsigned int> true_pathway(n_genes);
    unordered_map<string, unsigned int> counts;
    for (size_t i = 0; i < n_mc_samples; i++) {
        sample_pathway_uniform(random, n_pathways, true_pathway);
        string str = "";
        for (size_t g = 0; g < n_genes; g++) {
            if (g < n_genes - 1)
                str += to_string(true_pathway[g]) + "_";
            else
                str += to_string(true_pathway[g]);
        }
        counts[str] += 1;
    }
    
    double est = 0.0;
    for (auto it = counts.begin(); it != counts.end(); ++it) {
        est = (double)it->second/n_mc_samples;
        cout << it->first << ": " << est << endl;
        assert(abs(est - 1./6) < 0.01);
    }
    cout << "Prior sampling test passed!" << endl;
}

void test_log_marginal_estimates()
{
    // for small number of genes and pathway, we can compute the exact value of the marginal likelihood by
    // enumerating over all possibilities
    // then, we can check that the marginal likelihood estimate from IS/SMC approaches this value
    double fbp = 0.127;
    double bgp = 0.1721;
    double log_lik;

    unsigned int n_patients = 10;
    unsigned int n_pathways = 2;
    unsigned int n_genes = 3;

    unsigned int n_is_iter = 1;
    unsigned int n_is_particles = 10000;

    unsigned int n_smc_iter = 3;
    unsigned int n_smc_particles = 1000;
    double *exact_log_marginal_lik = new double[n_smc_iter-1];
    for (size_t i = 0; i < n_smc_iter; i++) {
        exact_log_marginal_lik[i] = DOUBLE_NEG_INF;
    }

    gsl_rng *random = generate_random_object(123);
    vector<unsigned int> true_pathway(n_genes);
    vector<unsigned int> stages(n_patients);
    
    unsigned int *pathway = new unsigned int[n_genes];
    
    sample_pathway_uniform(random, n_pathways, true_pathway);
    sample_stages_uniform(random, n_pathways, stages);
    gsl_matrix *data_matrix = simulate_data(random, n_pathways, true_pathway, stages);
    cout << "Clean matrix: " << endl;
    cout << "---" << endl;
    for (size_t i = 0; i < n_patients; i++) {
        for (size_t j = 0; j < n_genes; j++) {
            cout << gsl_matrix_get(data_matrix, i, j) << " ";
        }
        cout << endl;
    }
    cout << "---" << endl;
    cout << "Data matrix: " << endl;
    cout << "---" << endl;
    add_noise(random, fbp, bgp, data_matrix);
    for (size_t i = 0; i < n_patients; i++) {
        for (size_t j = 0; j < n_genes; j++) {
            cout << gsl_matrix_get(data_matrix, i, j) << " ";
        }
        cout << endl;
    }
    cout << "---" << endl;
    
    // enumerate over possible pathway using for loops (one per gene)
    unordered_map<string, double> state_log_liks;
    string pathway_str;
    for (size_t i = 0; i < n_pathways; i++) {
        pathway[0] = i;
        for (size_t j = 0; j < n_pathways; j++) {
            pathway[1] = j;
            for (size_t k = 0; k < n_pathways; k++) {
                pathway[2] = k;

                log_lik = compute_likelihood_from_matrix(data_matrix, pathway, n_pathways, n_genes, false, fbp, bgp);
                if (log_lik != DOUBLE_NEG_INF) {
                    pathway_str = to_string(pathway[0]) + ", " + to_string(pathway[1]) + ", " + to_string(pathway[2]);
                    cout << pathway_str << ", " << log_lik << endl;
                    state_log_liks[pathway_str] = log_lik;
                    for (size_t t = 1; t < n_smc_iter; t++) {
                        double beta_t = ((double)(t) / (n_smc_iter-1));
                        exact_log_marginal_lik[t-1] = log_add(exact_log_marginal_lik[t-1], beta_t * log_lik);
                    }
                }
            }
        }
    }

    for (auto it = state_log_liks.begin(); it != state_log_liks.end(); ++it) {
        state_log_liks[it->first] = exp(state_log_liks[it->first] - exact_log_marginal_lik[n_smc_iter-2]);
        cout << it->first << ", " << state_log_liks[it->first] << endl;
    }

    // add the prior factor
    double log_prior = log(1./6);
    for (size_t i = 0; i < n_smc_iter-1; i++) {
        exact_log_marginal_lik[i] += log_prior;
        if (i > 0) {
            cout << (i+1) << "/" << (n_smc_iter-1) << ", " << exact_log_marginal_lik[i] << ", " << (exact_log_marginal_lik[i] - exact_log_marginal_lik[i-1]) << endl;
        } else {
            cout << (i+1) << "/" << (n_smc_iter-1) << ", " << exact_log_marginal_lik[i] << ", " << exact_log_marginal_lik[i] << endl;
        }
    }

    long seed = 21;
    test_marginal_likelihood_estimate_smc(seed, data_matrix, n_pathways, n_is_particles, n_is_iter, fbp, bgp, state_log_liks, exact_log_marginal_lik[n_smc_iter-2]);
    test_marginal_likelihood_estimate_smc(seed, data_matrix, n_pathways, n_smc_particles, n_smc_iter, fbp, bgp, state_log_liks, exact_log_marginal_lik[n_smc_iter-2]);
}

int main()
{
    boost::filesystem::path curr_path = boost::filesystem::current_path();
    cout << "Exec dir: " <<  curr_path.string() << endl;
    //string data_path = "../../data/test";
    string data_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/test";

    test_likelihood(data_path);
    test_log_marginal_estimates();
    test_prior();
    return 0;
}

