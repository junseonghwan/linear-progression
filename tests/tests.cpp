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
#include <unordered_map>

#include <boost/filesystem.hpp>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "lpm.hpp"
#include "lpm_likelihood.hpp"

using namespace std;

// 1. test the likelihood calculation
void test_likelihood(string data_path, bool has_passenger)
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
    
    unsigned int idx;

    // free these objects within the loop
    gsl_matrix *data_matrix = 0;
    gsl_matrix *gold = 0;
    unsigned int *true_pathway = 0;

    for (unsigned int rep = 0; rep < 10; rep++) {
        rep_path = data_path + "/rep" + to_string(rep);
        gold_standard_file_name = rep_path + "/gold_standard.csv";
        data_file_name = rep_path + "/matrix.csv";
        pathway_file_name = rep_path + "/generative_mem_mat.csv";
        param_file_name = rep_path + "/parameters.csv";

        // read the data matrix
        data_matrix = read_csv(data_file_name);
        n_genes = data_matrix->size2;

        // read the gold standard likelihood file
        gold = read_csv(gold_standard_file_name, false);

        // read the true pathway
        true_pathway = new unsigned int[n_genes];
        n_pathways = read_ground_truth_pathway_from_matrix(pathway_file_name, true_pathway);
        if (has_passenger) {
            n_pathways--;
        }
        
        // read the parameters
        read_error_params(param_file_name, fbp, bgp);

        for (idx = 0; idx < gold->size1; idx++) {
            n_patients = gsl_matrix_get(gold, idx, 0);
            gold_log_lik = gsl_matrix_get(gold, idx, 1);

            gsl_matrix_view sub_matrix = gsl_matrix_submatrix(data_matrix, 0, 0, n_patients, n_genes);
            log_lik = compute_likelihood_from_matrix(&sub_matrix.matrix, true_pathway, n_pathways, n_genes, has_passenger, fbp, bgp);
            assert(abs(log_lik - gold_log_lik) < LIKELIHOOD_ERR_TOL);
            cout << n_patients << ": " << gold_log_lik << ", " << log_lik << endl;
        }
        
        // free the allocated memory
        delete [] true_pathway;
        delete data_matrix;
        delete gold;
    }
    
    cout << "Likelihood test passed!" << endl;
}

// run SMC to estimate the log marginal likelihood
//double test_marginal_likelihood_estimate_smc(long seed, gsl_matrix *data_matrix, unsigned int n_pathways, unsigned int n_particles, unsigned int n_smc_iter, unsigned int n_mcmc_iter, bool is_lik_tempered, double fbp, double bgp, unordered_map<string, double> true_probs, double exact_log_marginal_lik)
//{
//    double logZ;
//    double swap_prob = 0;
//    unsigned int n_genes = data_matrix->size2;
//
//    unsigned int *states = new unsigned int[n_particles * n_genes];
//    double *log_weights = new double[n_particles];
//    // run importance sampling with proposal from the prior
//    logZ = run_smc_from_matrix(seed, data_matrix, n_pathways, n_particles, n_smc_iter, n_mcmc_iter, false, swap_prob, fbp, bgp, states, log_weights, is_lik_tempered);
//    cout << "Estimated logZ: " << logZ << endl;
//    //assert(abs(exact_log_marginal_lik - logZ)/abs(exact_log_marginal_lik) < SMC_MARGINAL_ERR_TOL);
//
//    // perform tests using states and log_weights
//    string pathway_str;
//    unsigned int idx = 0;
//    unordered_map<string, double> counts;
//    normalize_destructively(log_weights, n_particles);
//    for (unsigned int i = 0; i < n_particles; i++) {
//        stringstream ss;
//        for (unsigned int j = 0; j < n_genes; j++) {
//            idx = i * n_genes + j;
//            ss << states[idx];
//            if (j < n_genes - 1) {
//                ss << ", ";
//            }
//        }
//        pathway_str = ss.str();
//        if (!counts.count(pathway_str)) {
//            counts[pathway_str] = 0;
//        }
//        counts[pathway_str] += log_weights[i];
//        //cout << pathway_str << ", " << counts[pathway_str] << endl;
//    }
//
//    double prob_est;
//    double prob_truth;
//    for (auto it = true_probs.begin(); it != true_probs.end(); ++it) {
//        pathway_str = it->first;
//        prob_est = counts[pathway_str];
//        prob_truth = true_probs[pathway_str];
//        cout << pathway_str << ", " << prob_est << ", " << prob_truth << endl;
//        assert(abs(prob_est - prob_truth) < SMC_MARGINAL_ERR_TOL);
//    }
//
//    cout << "SMC test passed!" << endl;
//    return logZ;
//
//    delete [] states;
//    delete [] log_weights;
//}

// for small number of genes and pathway, we can compute the exact value of the marginal likelihood by
// enumerating over all possibilities
// then, we can check that the marginal likelihood estimate from IS/SMC approaches this value
//void test_log_marginal_estimates()
//{
//
//    double fbp = 0.127;
//    double bgp = 0.1721;
//    double log_lik;
//
//    unsigned int n_patients = 10;
//    unsigned int n_pathways = 2;
//    unsigned int n_genes = 3;
//
//    unsigned int n_is_iter = 1;
//    unsigned int n_is_particles = 10000;
//
//    unsigned int n_smc_iter = 10;
//    unsigned int n_smc_particles = 20000;
//    double *exact_log_marginal_lik = new double[n_smc_iter];
//    for (unsigned int i = 0; i < n_smc_iter; i++) {
//        exact_log_marginal_lik[i] = DOUBLE_NEG_INF;
//    }
//
//    gsl_rng *random = generate_random_object(123);
//    vector<unsigned int> true_pathway(n_genes);
//    vector<unsigned int> stages(n_patients);
//
//    unsigned int *pathway = new unsigned int[n_genes];
//
//    propose_pathway(random, n_pathways, true_pathway);
//    sample_stages_uniform(random, n_pathways, stages);
//    gsl_matrix *data_matrix = simulate_data(random, n_pathways, true_pathway, stages);
//    add_noise(random, fbp, bgp, data_matrix);
//
//    // enumerate over possible pathway using for loops (one per gene)
//    unordered_map<string, double> state_log_liks;
//    string pathway_str;
//    for (unsigned int i = 0; i < n_pathways; i++) {
//        pathway[0] = i;
//        for (unsigned int j = 0; j < n_pathways; j++) {
//            pathway[1] = j;
//            for (unsigned int k = 0; k < n_pathways; k++) {
//                pathway[2] = k;
//
//                log_lik = compute_likelihood_from_matrix(data_matrix, pathway, n_pathways, n_genes, false, fbp, bgp);
//                if (log_lik != DOUBLE_NEG_INF) {
//                    pathway_str = to_string(pathway[0]) + ", " + to_string(pathway[1]) + ", " + to_string(pathway[2]);
//                    cout << pathway_str << ", " << log_lik << endl;
//                    state_log_liks[pathway_str] = log_lik;
//                    for (unsigned int t = 0; t < n_smc_iter; t++) {
//                        double beta_t = ((double)(t+1) / n_smc_iter);
//                        exact_log_marginal_lik[t] = log_add(exact_log_marginal_lik[t], beta_t * log_lik);
//                    }
//                }
//            }
//        }
//    }
//
//    for (auto it = state_log_liks.begin(); it != state_log_liks.end(); ++it) {
//        state_log_liks[it->first] = exp(state_log_liks[it->first] - exact_log_marginal_lik[n_smc_iter-1]);
//        cout << it->first << ", " << state_log_liks[it->first] << endl;
//    }
//
//    // add the prior factor
//    double log_prior = log(1./6);
//    for (unsigned int i = 0; i < n_smc_iter; i++) {
//        exact_log_marginal_lik[i] += log_prior;
//        if (i > 0) {
//            cout << (i+1) << "/" << n_smc_iter << ", " << exact_log_marginal_lik[i] << ", " << (exact_log_marginal_lik[i] - exact_log_marginal_lik[i-1]) << endl;
//        } else {
//            cout << (i+1) << "/" << n_smc_iter << ", " << exact_log_marginal_lik[i] << ", " << exact_log_marginal_lik[i] << endl;
//        }
//    }
//
//    long seed = gsl_rng_get(random);
//    test_marginal_likelihood_estimate_smc(seed, data_matrix, n_pathways, n_is_particles, n_is_iter, 1, false, fbp, bgp, state_log_liks, exact_log_marginal_lik[n_smc_iter-1]);
//
//    // test SMC algorithm using MH + tempering: n_mcmc_iter > 0 + is_lik_tempered = true
//    unsigned int n_mcmc_iter = 5;
//    unsigned int n_reps = 20;
//    bool is_lik_tempered = true;
//    vector<double> log_marginals(n_reps);
//    for (size_t i = 0; i < n_reps; i++) {
//        seed = gsl_rng_get(random);
//        double logZ = test_marginal_likelihood_estimate_smc(seed, data_matrix, n_pathways, n_smc_particles, n_smc_iter, n_mcmc_iter, is_lik_tempered, fbp, bgp, state_log_liks, exact_log_marginal_lik[n_smc_iter-1]);
//        log_marginals[i] = logZ;
//    }
//
//    double mean = gsl_stats_mean(log_marginals.data(), 1, log_marginals.size());
//    double sd = gsl_stats_sd_m(log_marginals.data(), 1, log_marginals.size(), mean);
//
//    // test that the 95% CI covers the truth
//    double lb = mean - 1.96 * sd;
//    double ub = mean + 1.96 * sd;
//    assert(lb <= exact_log_marginal_lik[n_smc_iter-1] && exact_log_marginal_lik[n_smc_iter-1] <= ub);
//    cout << "Truth: " << exact_log_marginal_lik[n_smc_iter - 1] << endl;
//    cout << "95% CI: " << "[" << lb << ", " << ub << "]" << endl;
//
//    delete [] exact_log_marginal_lik;
//    delete [] pathway;
//    delete random;
//}

void test_uniform_prior_calculation()
{
    double log_uniform_prior1 = -log_pathway_uniform_prior(2, 3);
    double log_uniform_prior2 = -log_pathway_uniform_prior_log_scale(2, 3);
    assert(abs(log_uniform_prior1 - log_uniform_prior2) < ERR_TOL);
    assert(abs(log_uniform_prior2 - log(1./6)) < ERR_TOL);

    log_uniform_prior1 = -log_pathway_uniform_prior(5, 50);
    log_uniform_prior2 = -log_pathway_uniform_prior_log_scale(5, 50);
    assert(abs(log_uniform_prior1 - log_uniform_prior2) < ERR_TOL);

    unsigned int n_genes = 3;
    unsigned int n_pathways = 2;
    double ret = -log_pathway_uniform_prior(n_pathways, n_genes);
    assert(ret == -log(6));

    n_genes = 8;
    n_pathways = 3;
    ret = -log_pathway_uniform_prior(n_pathways, n_genes);

    // sample pathways from the prior
    gsl_rng *random = generate_random_object(1);
    vector<unsigned int> pathway(n_genes);
    unordered_map<string, unsigned int> counts;
    unordered_map<string, double> probs;
    unsigned int n_mc_samples = 1000000;
    for (unsigned int n = 0; n < n_mc_samples; n++) {
        propose_pathway(random, n_pathways, pathway);
        string s = "";
        for (unsigned int g = 0; g < n_genes; g++) {
            s += to_string(pathway[g]);
            if (g < n_genes - 1) {
                s += ", ";
            }
        }
        if (!counts.count(s)) {
            counts[s] = 0;
            probs[s] = log_pathway_proposal(pathway, n_pathways);
        }
        counts[s] += 1;
    }

    cout << "Size: " << counts.size() << endl;
    assert(abs(ret - log(counts.size())));

    double prob = 0.0;
    for (auto it = counts.begin(); it != counts.end(); ++it) {
        prob = (double)it->second/n_mc_samples;
        ret = exp(probs[it->first]);
        cout << it->first << ", " << ret << ", " << prob << endl;
        assert(abs(prob - ret) < ERR_TOL);
    }
    
    delete random;
}

void test_run()
{
//    double error = 0.05;
//    const char* input_data = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment2/With_passengers/5/error0.05/rep0/matrix.csv";
//    const char* true_pathway_file = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment2/With_passengers/5/error0.05/rep0/generative_mem_mat.csv";
//    const char* output_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment2/With_passengers/5/error0.05/rep0/mcmc/";

    double error = 0.1;
    const char* input_data = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment1/With_passengers/Uniform/error0.1/rep0/matrix.csv";
    const char* true_pathway_file = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment1/With_passengers/Uniform/error0.1/rep0/generative_mem_mat.csv";
    const char* output_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment1/With_passengers/Uniform/error0.1/rep0/mcmc/";

    unsigned int *true_pathway = new unsigned int[100];
    read_ground_truth_pathway_from_matrix(true_pathway_file, true_pathway);
    double log_lik_at_truth = compute_likelihood(input_data, true_pathway, 5, 100, true, error, error);
    cout << log_lik_at_truth << endl;

    //double prior_passenger_prob = 0.95;
    //run_mcmc(121, input_data, output_path, 4, 20000, 20, 100, true, 0.1, 0.0, 0.2, 0.05, prior_passenger_prob);
    perform_prediction(output_path, input_data, 5, true, true_pathway);

//    gsl_rng *random = generate_random_object(121);
//    vector<double> log_Zs;
//    for (size_t i = 0; i < 5; i++) {
//        double log_Z = run_smc(gsl_rng_get(random), input_data, 6, 10000, 100, 1, true, 0.1, 0.05, 0.05, 0, 0, true, 6);
//        log_Zs.push_back(log_Z);
//        cout << log_Z << endl;
//    }
//    double log_Z_mean = gsl_stats_mean(log_Zs.data(), 1, log_Zs.size());
//    double log_Z_sd = gsl_stats_sd_m(log_Zs.data(), 1, log_Zs.size(), log_Z_mean);
//    double lb = log_Z_mean - 1.96 * log_Z_sd;
//    double ub = log_Z_mean + 1.96 * log_Z_sd;
//    cout << "95% CI: [" << lb << ", " << ub << "]" <<  endl;
//    gsl_rng_free(random);
    
    //model_selection(1, input_data, 5, 1, 500, 100, 1, true, 0.1, 0, 0, 1, 6, 0, 0, true);
}

int main(int argc, char *argv[])
{
    bool run_test = false;
    if (run_test) {
        if (argc != 2)
            exit(-1);
        
        string root_path(argv[1]);
        string lik_test_path1 = root_path + "/data/test/Without_passengers/";
        string lik_test_path2 = root_path + "/data/test/With_passengers/";

        test_uniform_prior_calculation();
        test_likelihood(lik_test_path1, false);
        test_likelihood(lik_test_path2, true);
        //test_log_marginal_estimates();
    } else {
        test_run();
    }

    return 0;
}

