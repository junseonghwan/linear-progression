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
    unsigned int idx = 0;
    unordered_map<string, double> counts;
    normalize_destructively(log_weights, n_particles);
    for (unsigned int i = 0; i < n_particles; i++) {
        stringstream ss;
        for (unsigned int j = 0; j < n_genes; j++) {
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
        //cout << pathway_str << ", " << counts[pathway_str] << endl;
    }

    double prob_est;
    double prob_truth;
    for (auto it = true_probs.begin(); it != true_probs.end(); ++it) {
        pathway_str = it->first;
        prob_est = counts[pathway_str];
        prob_truth = true_probs[pathway_str];
        //cout << pathway_str << ", " << prob_est << ", " << prob_truth << endl;
        assert(abs(prob_est - prob_truth) < SMC_MARGINAL_ERR_TOL);
    }
    
    cout << "SMC test passed!" << endl;
    
    delete [] states;
    delete [] log_weights;
}

// for small number of genes and pathway, we can compute the exact value of the marginal likelihood by
// enumerating over all possibilities
// then, we can check that the marginal likelihood estimate from IS/SMC approaches this value
void test_log_marginal_estimates()
{

    double fbp = 0.127;
    double bgp = 0.1721;
    double log_lik;

    unsigned int n_patients = 10;
    unsigned int n_pathways = 2;
    unsigned int n_genes = 3;

    unsigned int n_is_iter = 1;
    unsigned int n_is_particles = 10000;

    unsigned int n_smc_iter = 3;
    unsigned int n_smc_particles = 100000;
    double *exact_log_marginal_lik = new double[n_smc_iter-1];
    for (unsigned int i = 0; i < n_smc_iter; i++) {
        exact_log_marginal_lik[i] = DOUBLE_NEG_INF;
    }

    gsl_rng *random = generate_random_object(123);
    vector<unsigned int> true_pathway(n_genes);
    vector<unsigned int> stages(n_patients);
    
    unsigned int *pathway = new unsigned int[n_genes];
    
    sample_pathway_from_prior(random, n_pathways, true_pathway);
    sample_stages_uniform(random, n_pathways, stages);
    gsl_matrix *data_matrix = simulate_data(random, n_pathways, true_pathway, stages);
    cout << "Clean matrix: " << endl;
    cout << "---" << endl;
    for (unsigned int i = 0; i < n_patients; i++) {
        for (unsigned int j = 0; j < n_genes; j++) {
            cout << gsl_matrix_get(data_matrix, i, j) << " ";
        }
        cout << endl;
    }
    cout << "---" << endl;
    cout << "Data matrix: " << endl;
    cout << "---" << endl;
    add_noise(random, fbp, bgp, data_matrix);
    for (unsigned int i = 0; i < n_patients; i++) {
        for (unsigned int j = 0; j < n_genes; j++) {
            cout << gsl_matrix_get(data_matrix, i, j) << " ";
        }
        cout << endl;
    }
    cout << "---" << endl;
    
    // enumerate over possible pathway using for loops (one per gene)
    unordered_map<string, double> state_log_liks;
    string pathway_str;
    for (unsigned int i = 0; i < n_pathways; i++) {
        pathway[0] = i;
        for (unsigned int j = 0; j < n_pathways; j++) {
            pathway[1] = j;
            for (unsigned int k = 0; k < n_pathways; k++) {
                pathway[2] = k;

                log_lik = compute_likelihood_from_matrix(data_matrix, pathway, n_pathways, n_genes, false, fbp, bgp);
                if (log_lik != DOUBLE_NEG_INF) {
                    pathway_str = to_string(pathway[0]) + ", " + to_string(pathway[1]) + ", " + to_string(pathway[2]);
                    cout << pathway_str << ", " << log_lik << endl;
                    state_log_liks[pathway_str] = log_lik;
                    for (unsigned int t = 1; t < n_smc_iter; t++) {
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
    for (unsigned int i = 0; i < n_smc_iter-1; i++) {
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
    
    delete [] exact_log_marginal_lik;
    delete [] pathway;
    delete random;
}

void test_prior_sampling()
{
    gsl_rng *random = generate_random_object(123);

    unsigned int *genes = new unsigned int[4];
    unsigned int *dest = new unsigned int[3];
    for (unsigned int i = 0; i < 4; i++) {
        genes[i] = i;
    }
    unordered_map<string, unsigned int> map;
    for (unsigned int i = 0; i < 100000; i++) {
        gsl_ran_choose(random, dest, 3, genes, 4, sizeof(unsigned int));
        string str = to_string(dest[0]) + ", " + to_string(dest[1]) + ", " + to_string(dest[2]);
        if (map.count(str) == 0) {
            map[str] = 0;
        }
        map[str] += 1;
    }
    for (auto it = map.begin(); it != map.end(); ++it) {
        cout << it->first << ": " << (double)map[it->first]/100000 << endl;
    }

    unsigned int n_pathways = 2;
    unsigned int n_genes = 3;
    unsigned int n_mc_samples = 100000;
    
    vector<unsigned int> true_pathway(n_genes);
    unordered_map<string, unsigned int> counts;
    for (unsigned int i = 0; i < n_mc_samples; i++) {
        sample_pathway_from_prior(random, n_pathways, true_pathway);
        string str = "";
        for (unsigned int g = 0; g < n_genes; g++) {
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
    
    delete random;
    delete [] genes;
    delete [] dest;
}

void test_prior_calculation()
{
    unsigned int n_genes = 3;
    unsigned int n_pathways = 2;
    double ret = log_pathway_uniform_prior(n_pathways, n_genes);
    assert(ret == -log(6));

    n_genes = 4;
    n_pathways = 2;
    ret = log_pathway_uniform_prior(n_pathways, n_genes);

    // sample pathways from the prior
    gsl_rng *random = generate_random_object(1);
    vector<unsigned int> pathway(n_genes);
    unordered_map<string, unsigned int> counts;
    unordered_map<string, double> probs;
    unsigned int n_mc_samples = 1000000;
    for (unsigned int n = 0; n < n_mc_samples; n++) {
        sample_pathway_from_prior(random, n_pathways, pathway);
        string s = "";
        for (unsigned int g = 0; g < n_genes; g++) {
            s += to_string(pathway[g]);
            if (g < n_genes - 1) {
                s += ", ";
            }
        }
        if (!counts.count(s)) {
            counts[s] = 0;
            probs[s] = log_pathway_prior(pathway, n_pathways);
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

void test_posterior()
{
    gsl_rng *random = generate_random_object(1);

    long seed;
    double error_max = 0.2;
    unsigned int num_reps = 20;
    unsigned int n_pathways = 3;
    unsigned int n_patients = 200;
    unsigned int n_genes = 5;

    unsigned int n_pg_iter = 50;
    unsigned int n_particles = 200;
    unsigned int n_smc_iter = 2*n_genes;
    unsigned int n_kernel_iter = 1;
    unsigned int n_mh_w_gibbs_iter = 20;
    bool has_passenger = false;
    double swap_prob = 0.2;

    vector<unsigned int> stages(n_patients);
    vector<unsigned int> true_pathway(n_genes);
    vector<unsigned int> row_sum(n_patients);

    sample_stages_uniform(random, n_pathways, stages);
    sample_pathway_from_prior(random, n_pathways, true_pathway);
    
    double bgp, fbp, chi2 = 0.0;
    for (unsigned int rep = 0; rep < num_reps; rep++) {
        // 1. sample parameters
        // 2. add noise to the data
        // 3. run particle Gibbs to generate posterior over the parameters
        // 4. compute quantile of the sampled parametrs amongst the posterior samples: q_i
        // 5. convert the quantiles to standard normal: Z_i = \Phi(q_i)
        bgp = gsl_ran_flat(random, 0.01, error_max);
        //fbp = gsl_ran_flat(random, 0.01, error_max);;
        fbp = bgp;
        gsl_matrix *data_matrix = simulate_data(random, n_pathways, true_pathway, stages);
        add_noise(random, fbp, bgp, data_matrix);
        
        // compute the likelihood at the true parameters
        compute_row_sum(*data_matrix, row_sum);
        LinearProgressionParameters params(fbp, bgp);
        //LinearProgressionState true_state(n_genes, n_pathways, false);
        LinearProgressionState true_state(*data_matrix, row_sum, n_genes, n_pathways, false);
        for (unsigned int idx = 0; idx < n_genes; idx++) {
            true_state.update_pathway_membership(idx, true_pathway[idx]);
        }
        double true_log_lik = compute_pathway_likelihood(true_state, params);
        cout << "True log lik: " << true_log_lik << endl;

        seed = gsl_rng_get(random);
        ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg = run_pg_from_matrix(seed, data_matrix, n_pathways, n_pg_iter, n_particles, n_smc_iter, n_kernel_iter, n_mh_w_gibbs_iter, has_passenger, swap_prob, 0.0, error_max);

        //vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &ret_states = pg.get_states();
        //vector<double> log_marginals = pg.get_log_marginal_likelihoods();
        vector<shared_ptr<LinearProgressionParameters>> &ret_params = pg.get_parameters();
        
        unsigned int bgp_q = 0;
        for (unsigned int j = 0; j < n_pg_iter; j++) {
            LinearProgressionParameters &param = *ret_params[j].get();
            if (bgp > param.get_bgp()) {
                bgp_q += 1;
            }
        }
        double u =  (double)bgp_q/n_pg_iter;
        double zi = gsl_cdf_gaussian_Qinv(u, 1);
        chi2 += pow(zi, 2);
        cout << "BGP: " << bgp << endl;
        cout << "FBP: " << fbp << endl;
        cout << "Z_i: " << zi << endl;
        
        delete data_matrix;
    }
    // 6. test that the mean of Z = 0 using {Z_i}. If the null hypothesis is rejected, then the posterior code is incorrectly implemented.
    cout << "chi^2: " << chi2 << endl;
    // check that |z| < 1.96
    double pval = 1 - gsl_cdf_chisq_P(chi2, num_reps);
    assert(pval > 0.05);
}

int main()
{
    boost::filesystem::path curr_path = boost::filesystem::current_path();
    cout << "Exec dir: " <<  curr_path.string() << endl;
    //string data_path = "../../data/test";
    string data_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/test";

    test_likelihood(data_path);
    test_prior_sampling();
    test_prior_calculation();
    test_log_marginal_estimates();
    test_posterior(); // this takes about 5 minutes depending on the hardware
    return 0;
}

