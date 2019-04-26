//
//  lpm.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-28.
//

#include "lpm.hpp"

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>

#include <gsl/gsl_sf.h>

#include <omp.h>

#include "data_util.hpp"
#include "lpm_likelihood.hpp"
#include "lpm_model.hpp"

using namespace std;

void run_mcmc(long seed,
              const char *dat_file,
              const char *output_path,
              unsigned int *initial_state,
              unsigned int n_mcmc_iter,
              double swap_prob,
              double fbp_max,
              double bgp_max,
              double mh_proposal_sd)
{
    // run MH on the states and params
    
}

void run_pg(long seed,
            const char *dat_file,
            const char *output_path,
            unsigned int model_len,
            unsigned int n_pg_iter,
            unsigned int n_particles,
            unsigned int n_smc_iter,
            unsigned int n_kernel_iter,
            unsigned int n_mh_w_gibbs_iter,
            bool has_passenger,
            double swap_prob,
            double fbp_max,
            double bgp_max,
            double mh_proposal_sd)
{
    // convert char* to string
    string dat_file_str(dat_file, strlen(dat_file));

    // load the data
    gsl_matrix *obs_matrix = read_csv(dat_file_str, false);
    unsigned int n_genes = obs_matrix->size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    run_pg_from_matrix(seed, obs_matrix, model_len, n_pg_iter, n_particles, n_smc_iter, n_kernel_iter, n_mh_w_gibbs_iter, has_passenger, swap_prob, fbp_max, bgp_max, 0, 0, mh_proposal_sd, output_path);

    delete obs_matrix;
}

ParticleGibbs<LinearProgressionState, LinearProgressionParameters> run_pg_from_matrix(
                                                                                      long seed,
                                                                                      gsl_matrix *obs_matrix,
                                                                                      unsigned int model_len,
                                                                                      unsigned int n_pg_iter,
                                                                                      unsigned int n_particles,
                                                                                      unsigned int n_smc_iter,
                                                                                      unsigned int n_kernel_iter,
                                                                                      unsigned int n_mh_w_gibbs_iter,
                                                                                      bool has_passenger,
                                                                                      double swap_prob,
                                                                                      double fbp_max,
                                                                                      double bgp_max,
                                                                                      vector<double> *fbps,
                                                                                      vector<double> *bgps,
                                                                                      double mh_proposal_sd,
                                                                                      const char *output_path)
{
    unsigned int n_patients = obs_matrix->size1;
    unsigned int n_genes = obs_matrix->size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    // output the settings
    cout << "Running PG {" << endl;
    cout << "\tNum patients: " << n_patients << endl;
    cout << "\tModel length: " << model_len << endl;
    cout << "\tNum PG iter: " << n_pg_iter << endl;
    cout << "\tNum particles: " << n_particles << endl;
    cout << "\tNum SMC iter: " << n_smc_iter << endl;
    cout << "\tNum kernel iter: " << n_kernel_iter << endl;
    cout << "\tNum MHwGibbs iter: " << n_mh_w_gibbs_iter << endl;
    cout << "\tAllocate passenger pathway: " << to_string(has_passenger) << endl;
    cout << "\tSwap prob: " << swap_prob << endl;
    cout << "\tFBP max: " << fbp_max << endl;
    cout << "\tBGP max: " << bgp_max << endl;
    cout << "}" << endl;

    // compute the row sum for each patient
    vector<unsigned int> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);

    // allocate random object
    gsl_rng *random = generate_random_object(seed);

    // smc options
    SMCOptions smc_options;
    smc_options.num_particles = n_particles;
    smc_options.ess_threshold = 1.0;
    smc_options.resample_last_round = false;
    smc_options.main_seed = gsl_rng_get(random);
    smc_options.resampling_seed = gsl_rng_get(random);

    // pmcmc options
    PMCMCOptions pmcmc_options(gsl_rng_get(random), n_pg_iter);

    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, false);

    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);

    // Declare param proposal object
    LPMParamProposal pg_proposal(n_mh_w_gibbs_iter, fbp_max, bgp_max, mh_proposal_sd);

    // Initialize PG object
    ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(pmcmc_options, csmc, pg_proposal);
    
    // run (output timing in seconds)
    auto start = std::chrono::high_resolution_clock::now();
    pg.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds elapsed." <<  endl;

    if (output_path != 0) {
        // output the states and parameters
        string output_path_str(output_path, strlen(output_path));
        write_pg_output(output_path_str, pg, pg_proposal);
    }
    if (fbps != 0) {
        *fbps = pg_proposal.get_fbps();
    }
    if (bgps != 0) {
        *bgps = pg_proposal.get_bgps();
    }
    return pg;
    
    delete random;
}

double run_smc(long seed,
                const char *dat_file,
                unsigned int model_len,
                unsigned int n_particles,
                unsigned int n_smc_iter,
                unsigned int n_kernel_iter,
                bool has_passenger,
                double swap_prob,
                double fbp,
                double bgp,
                unsigned int *states,
                double *log_weights,
                bool is_lik_tempered)
{
    // convert char* to string
    string dat_file_str(dat_file, strlen(dat_file));
    
    // load the data
    gsl_matrix *obs_matrix = read_csv(dat_file_str, false);
    unsigned int n_patients = obs_matrix->size1;

    // output the settings
    cout << "Running SMC {" << endl;
    cout << "\tInput data: " << dat_file_str << endl;
    cout << "\tNum patients: " << n_patients << endl;
    cout << "\tModel length: " << model_len << endl;
    cout << "\tNum particles: " << n_particles << endl;
    cout << "\tNum SMC iter: " << n_smc_iter << endl;
    cout << "\tNum kernel iter: " << n_kernel_iter << endl;
    cout << "\tAllocate passenger pathway: " << to_string(has_passenger) << endl;
    cout << "\tSwap prob: " << swap_prob << endl;
    cout << "\tFBP: " << fbp << endl;
    cout << "\tBGP: " << bgp << endl;
    cout << "}" << endl;

    return run_smc_from_matrix(seed, obs_matrix, model_len, n_particles, n_smc_iter, n_kernel_iter, has_passenger, swap_prob, fbp, bgp, states, log_weights, is_lik_tempered);
    
    delete obs_matrix;
}

double run_smc_from_matrix(long seed,
                           gsl_matrix *obs_matrix,
                           unsigned int model_len,
                           unsigned int n_particles,
                           unsigned int n_smc_iter,
                           unsigned int n_kernel_iter,
                           bool has_passenger,
                           double swap_prob,
                           double fbp,
                           double bgp,
                           unsigned int *states,
                           double *log_weights,
                           bool is_lik_tempered)
{
    unsigned int n_patients = obs_matrix->size1;
    unsigned int n_genes = obs_matrix->size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    // compute the row sum for each patient
    vector<unsigned int> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);

    // allocate random object
    gsl_rng *random = generate_random_object(seed);

    // smc options
    SMCOptions smc_options;
    smc_options.num_particles = n_particles;
    smc_options.ess_threshold = 1.0;
    smc_options.resample_last_round = false;
    smc_options.main_seed = gsl_rng_get(random);
    smc_options.resampling_seed = gsl_rng_get(random);

    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, is_lik_tempered);

    // Construct LPMParam object
    LinearProgressionParameters param(fbp, bgp);

    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> smc(smc_model, smc_options);

    // run (output timing in seconds)
    auto start = std::chrono::high_resolution_clock::now();
    smc.initialize(param);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds elapsed." <<  endl;

    if (states != 0 && log_weights != 0) {
        // fill the states and log_weights
        // dimension of states should be n_particles x n_genes
        // dimension of the log_weights should be n_particles
        unsigned int idx = 0;
        for (unsigned int i = 0; i < n_particles; i++) {
            log_weights[i] = smc.get_log_weight(n_smc_iter - 1, i);
            const LinearProgressionState &state = smc.get_state(n_smc_iter - 1, i);
            const vector<unsigned int> pathway = state.get_pathway_membership();
            for (unsigned int j = 0; j < n_genes; j++) {
                idx = i*n_genes + j;
                states[idx] = pathway[j];
            }
        }
    }
    
//    for (size_t i = 0; i < n_smc_iter; i++) {
//        cout << "Z" << i << ": " << (smc.get_log_norm(i) - log(n_particles)) << endl;
//    }

    return smc.get_log_marginal_likelihood();
    
    delete random;
}

double model_selection(long seed,
                       const char *dat_file,
                       unsigned int model_len,
                       unsigned int n_mc_samples,
                       unsigned int n_particles,
                       unsigned int n_smc_iter,
                       unsigned int n_kernel_iter,
                       bool has_passenger,
                       double swap_prob,
                       const double *fbps,
                       const double *bgps,
                       unsigned int n_threads,
                       double *ret_sum,
                       double *ret_smc)
{
    // convert char* to string
    string dat_file_str(dat_file, strlen(dat_file));
    
    // load the data
    gsl_matrix *obs_matrix = read_csv(dat_file_str, false);
    unsigned int n_patients = obs_matrix->size1;
    unsigned int n_genes = obs_matrix->size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    // output the settings
    cout << "Running Model Selection {" << endl;
    cout << "\tInput data: " << dat_file_str << endl;
    cout << "\tNum patients: " << n_patients << endl;
    cout << "\tModel length: " << model_len << endl;
    cout << "\tNum samples for parameters: " << n_mc_samples << endl;
    cout << "\tNum particles: " << n_particles << endl;
    cout << "\tNum SMC iter: " << n_smc_iter << endl;
    cout << "\tNum kernel iter: " << n_kernel_iter << endl;
    cout << "\tAllocate passenger pathway: " << to_string(has_passenger) << endl;
    cout << "\tSwap prob: " << swap_prob << endl;
    cout << "}" << endl;
    
    // compute the row sum for each patient
    vector<unsigned int> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);
    
    // allocate random object
    gsl_rng *random = generate_random_object(seed);
    
    // generate seeds and parameters
    vector<long> main_seeds, resampling_seeds;
    for (unsigned int i = 0; i < n_mc_samples; i++) {
        main_seeds.push_back(gsl_rng_get(random));
        resampling_seeds.push_back(gsl_rng_get(random));
    }
    
    vector<double> log_marginals(n_mc_samples);
    vector<double> smc_log_marginals(n_mc_samples);
    auto start = std::chrono::high_resolution_clock::now();
    
    printf("n_unique_states, fbp, bgp, log_marginal_sum, log_marginal_smc\n");
    omp_set_num_threads(n_threads);
#pragma omp parallel for
    for (unsigned int n = 0; n < n_mc_samples; n++) {
        // smc options
        SMCOptions smc_options;
        smc_options.num_particles = n_particles;
        smc_options.ess_threshold = 1.0;
        smc_options.resample_last_round = false;
        smc_options.main_seed = main_seeds[n];
        smc_options.resampling_seed = resampling_seeds[n];

        // LPM model proposal
        LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, 0, *obs_matrix, row_sum, swap_prob, has_passenger, false);

        // Declare conditional SMC object
        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);

        // Construct LPMParam object
        LinearProgressionParameters param(fbps[n], bgps[n]);

        // run SMC
        csmc.initialize(param);
        
        // get log marginal estimate from SMC
        if (ret_smc != 0) {
            ret_smc[n] = csmc.get_log_marginal_likelihood();
        }

        // enumerate over all particles that was ever generated to approximate \sum_x p(y | x, \theta)
        unordered_set<string> unique_states;
        log_marginals[n] = DOUBLE_NEG_INF;
        double log_prior = 0.0;
        for (unsigned int r = 0; r < n_smc_iter; r++) {
            for (unsigned int i = 0; i < n_particles; i++) {
                const LinearProgressionState &state = csmc.get_state(r, i);
                string key = state.to_string();
                if (unique_states.count(key) == 0) {
                    log_prior = log_pathway_proposal(state, n_genes);
                    log_marginals[n] = log_add(log_marginals[n], state.get_log_lik() + log_prior);
                    unique_states.insert(key);
                }
            }
        }
        if (ret_sum != 0) {
            ret_sum[n] = log_marginals[n];
        }

        printf("%zd, %f, %f, %f, %f\n", unique_states.size(), fbps[n], bgps[n], log_marginals[n], ret_smc[n]);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds elapsed." <<  endl;

    // combine log_marginals to approximate:
    // \int_{\theta} \sum_x p(y|x, k, \theta) p(x|k) p(\theta) d(\theta)
    double log_f_hat = DOUBLE_NEG_INF;
    for (unsigned int n = 0; n < n_mc_samples; n++) {
        log_f_hat = log_add(log_f_hat, log_marginals[n]);
    }

    return log_f_hat;
    
    delete obs_matrix;
    delete random;
}

double compute_likelihood(const char *dat_file,
                          unsigned int *pathway,
                          unsigned int model_len,
                          unsigned int n_genes,
                          bool has_passenger,
                          double fbp,
                          double bgp)
{
    string dat_file_str(dat_file, strlen(dat_file));
    
    // load the data
    gsl_matrix *obs_matrix = read_csv(dat_file_str, false);
    return compute_likelihood_from_matrix(obs_matrix, pathway, model_len, n_genes, has_passenger, fbp, bgp);
    
    delete obs_matrix;
}

double compute_likelihood_from_matrix(gsl_matrix *obs_matrix,
                                      unsigned int *pathway,
                                      unsigned int n_drivers,
                                      unsigned int n_genes,
                                      bool has_passenger,
                                      double fbp,
                                      double bgp)
{
    if (obs_matrix->size2 != n_genes)
    {
        cerr << "Error: Number of genes specified is " << n_genes << ". But the data contains " << obs_matrix->size2 << " genes." << endl;
        exit(-1);
    }
 
    unsigned int n_patients = obs_matrix->size1;

    // compute the row sum for each patient
    vector<unsigned int> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);

    // construct LPState from the pathway
    LinearProgressionState state(*obs_matrix, row_sum, n_genes, n_drivers, has_passenger);
    for (unsigned int i = 0; i < n_genes; i++) {
        state.update_pathway_membership(i, pathway[i]);
    }

    // construct params
    LinearProgressionParameters params(fbp, bgp);

    double log_lik = compute_pathway_likelihood(state, params);
    return log_lik;

    delete obs_matrix;
}

void generate_data(long seed,
                   const char *output_path,
                   unsigned int model_len,
                   unsigned int n_genes,
                   unsigned int n_patients,
                   double fbp,
                   double bgp)
{
    string output_path_str(output_path, strlen(output_path));

    gsl_rng *random = generate_random_object(seed);
    
    vector<unsigned int> stages(n_patients);
    vector<unsigned int> true_pathway(n_genes);
    
    sample_stages_uniform(random, model_len, stages);
    propose_pathway(random, model_len, true_pathway);
    gsl_matrix *data_matrix = simulate_data(random, model_len, true_pathway, stages);

    // output: stages, pathway, parameters, data matrix before contamination
    string stages_file = output_path_str + "/stages.csv";
    string pathway_file = output_path_str + "/true_pathways.csv";
    string clean_matrix_file = output_path_str + "/clean_matrix.csv";
    string data_matrix_file = output_path_str + "/data_matrix.csv";

    write_vector(stages_file, stages);
    write_vector(pathway_file, true_pathway);
    write_matrix_as_csv(clean_matrix_file, *data_matrix);
    
    add_noise(random, fbp, bgp, data_matrix);

    // output: data after contamination
    write_matrix_as_csv(data_matrix_file, *data_matrix);
    
    delete random;
    delete data_matrix;
}
