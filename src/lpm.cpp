//
//  lpm.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-28.
//

#include "lpm.hpp"

#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include "data_util.hpp"
#include "lpm_likelihood.hpp"
#include "lpm_model.hpp"

using namespace std;

bool do_swap_move(unsigned int num_driver_pathways, double unif, double pathway_swap_prob)
{
    return (num_driver_pathways >= 2 && unif < pathway_swap_prob);
}

void mh_move_params(gsl_rng *random,
                    LinearProgressionParameters &params,
                    LinearProgressionState &state,
                    unsigned int n_mh_w_gibbs_iter,
                    double mh_proposal_sd,
                    double error_max)
{
    double old_log_lik = compute_pathway_likelihood(state, params);
    double old_error_prob = params.get_bgp();
    for (size_t i = 0; i < n_mh_w_gibbs_iter; i++) {
        double new_error_prob = gsl_ran_gaussian(random, mh_proposal_sd) + old_error_prob;
        if (new_error_prob > 0 && new_error_prob <= error_max) {
            params.set_bgp(new_error_prob);
            params.set_fbp(new_error_prob);
            double new_log_lik = compute_pathway_likelihood(state, params);
            double log_u = log(gsl_ran_flat(random, 0.0, 1.0));
            if (log_u < (new_log_lik - old_log_lik)) {
                // accept
                old_log_lik = new_log_lik;
                old_error_prob = new_error_prob;
            } else {
                // revert
                params.set_bgp(old_error_prob);
                params.set_fbp(old_error_prob);
            }
        }
    }
}

unsigned int *pathway_indices = 0;
void mh_move_pathway(gsl_rng *random,
             LinearProgressionParameters &params,
             LinearProgressionState &state,
             unsigned int n_genes,
             double pathway_swap_prob,
             double prior_passenger_prob)
{
    double log_accept_ratio;
    double prev_log_lik = compute_pathway_likelihood(state, params);
    double new_log_lik;

    unsigned int *pathways_to_swap = new unsigned int[2];
    double log_unif = 0.0;
    bool is_swap_move = false;
    unsigned gene_idx, pathway1, pathway2;

    if (prev_log_lik == DOUBLE_NEG_INF) {
        cerr << "Error: previous log lik is neg inf! This means empty pathway was sampled in the previous step." << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < 3; i++) {
        double unif = gsl_ran_flat(random, 0.0, 1.0);
        is_swap_move = do_swap_move(state.get_num_driver_pathways(), unif, pathway_swap_prob);
        if (is_swap_move) {
            gsl_ran_choose(random, pathways_to_swap, 2, pathway_indices, state.get_num_driver_pathways(), sizeof(unsigned int));
            pathway1 = pathways_to_swap[0];
            pathway2 = pathways_to_swap[1];
            state.swap_pathways(pathway1, pathway2);
        } else {
            // sample a gene at random and sample pathway at random
            gene_idx = gsl_rng_uniform_int(random, n_genes);
            pathway1 = state.get_pathway_membership_of(gene_idx);
            pathway2 = gsl_rng_uniform_int(random, state.get_num_pathways());
            if (pathway1 == pathway2)
                continue;
            state.update_pathway_membership(gene_idx, pathway2);
        }

        new_log_lik = compute_pathway_likelihood(state, params);

        log_unif = log(gsl_ran_flat(random, 0.0, 1.0));
        log_accept_ratio = (new_log_lik - prev_log_lik);

        if (!is_swap_move) {
            if (pathway1 < state.get_num_driver_pathways() && pathway2 == state.get_num_driver_pathways()) {
                // assigning a driver gene to a passenger gene
                log_accept_ratio += (log(prior_passenger_prob) - log(1-prior_passenger_prob));
            } else if (pathway1 == state.get_num_driver_pathways() && pathway2 < state.get_num_driver_pathways()) {
                // assigning from passenger to driver
                log_accept_ratio += (log(1-prior_passenger_prob) - log(prior_passenger_prob));
            }
        }

        if (log_unif < log_accept_ratio) {
            // accept
            state.set_log_lik(new_log_lik);
            prev_log_lik = new_log_lik;
        } else {
            // revert
            if (is_swap_move) {
                state.swap_pathways(pathway2, pathway1);
            } else {
                state.update_pathway_membership(gene_idx, pathway1);
            }
            state.set_log_lik(prev_log_lik);
        }
    }
}

void run_mcmc(long seed,
              const char *dat_file,
              const char *output_path,
              unsigned int model_len,
              unsigned int n_mcmc_iter,
              unsigned int n_mh_w_gibbs_iter,
              unsigned int thinning,
              bool has_passenger,
              double swap_prob,
              double fbp_max,
              double bgp_max,
              double mh_proposal_sd,
              double prior_passenger_prob)
{
    // convert char* to string
    string dat_file_str(dat_file, strlen(dat_file));
    
    // load the data
    gsl_matrix *obs_matrix = read_csv(dat_file_str, false);
    unsigned int n_patients = obs_matrix->size1;
    unsigned int n_genes = obs_matrix->size2;

    // output the settings
    cout << "Running MCMC sampler {" << endl;
    cout << "\tInput data: " << dat_file_str << endl;
    cout << "\tNum patients: " << n_patients << endl;
    cout << "\tModel length: " << model_len << endl;
    cout << "\tNum MCMC iterations: " << n_mcmc_iter << endl;
    cout << "\tNum MHwGibbs iter: " << n_mh_w_gibbs_iter << endl;
    cout << "\tAllocate passenger pathway: " << to_string(has_passenger) << endl;
    cout << "\tSwap prob: " << swap_prob << endl;
    cout << "\tPrior passenger prob: " << prior_passenger_prob << endl;
    cout << "\tFBP max: " << fbp_max << endl;
    cout << "\tBGP max: " << bgp_max << endl;
    cout << "}" << endl;

    size_t burn_in = n_mcmc_iter / 10;
    size_t n_samples = (n_mcmc_iter-burn_in)/thinning;
    
    gsl_matrix *states = gsl_matrix_alloc(n_samples, n_genes);
    vector<double> bgps;
    vector<double> fbps;

    run_mcmc_from_matrix(seed, *obs_matrix, model_len, n_mcmc_iter, n_mh_w_gibbs_iter, thinning, burn_in, has_passenger, swap_prob, fbp_max, bgp_max, mh_proposal_sd, states, bgps, fbps, prior_passenger_prob);

    // print the samples to files
    string output_path_str(output_path, strlen(output_path));
    write_vector(output_path_str + "/bgps.csv" , bgps);
    write_vector(output_path_str + "/fbps.csv" , fbps);
    write_matrix_as_csv(output_path_str + "/states.csv" , *states);

    process_samples(output_path_str, *obs_matrix, states, bgps, fbps, model_len, has_passenger, prior_passenger_prob);
    
    gsl_matrix_free(states);

}

void run_mcmc_from_matrix(long seed,
                          const gsl_matrix &obs_matrix,
                          unsigned int model_len,
                          unsigned int n_mcmc_iter,
                          unsigned int n_mh_w_gibbs_iter,
                          unsigned int thinning,
                          unsigned int burn_in,
                          bool has_passenger,
                          double swap_prob,
                          double fbp_max,
                          double bgp_max,
                          double mh_proposal_sd,
                          gsl_matrix *states,
                          vector<double> &bgps,
                          vector<double> &fbps,
                          double prior_passenger_prob)
{
    unsigned int n_patients = obs_matrix.size1;
    unsigned int n_genes = obs_matrix.size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    vector<unsigned int> row_sum(n_patients);
    compute_row_sum(obs_matrix, row_sum);

    // allocate random object
    gsl_rng *random = generate_random_object(seed);

    pathway_indices = new unsigned int[model_len];
    for (unsigned int i = 0; i < model_len; i++) {
        pathway_indices[i] = i;
    }

    double bgp_init = gsl_ran_flat(random, 0, bgp_max);
    double fbp_init = bgp_init;
    if (fbp_max > 0) {
        fbp_init = gsl_ran_flat(random, 0, fbp_max);
    }
    double path_swap_prob = 0.1;
    LinearProgressionParameters params(fbp_init, bgp_init);
    LinearProgressionState *state = new LinearProgressionState(obs_matrix, row_sum, n_genes, model_len, has_passenger);
    //state->sample_min_valid_pathway(random);
    state->sample_pathway(random);
    state->set_log_lik(compute_pathway_likelihood(*state, params));
    double best_log_lik = DOUBLE_NEG_INF;

    size_t sample_idx = 0;
    for (size_t i = 0; i < n_mcmc_iter; i++) {
        mh_move_pathway(random, params, *state, n_genes, path_swap_prob, prior_passenger_prob);
        mh_move_params(random, params, *state, n_mh_w_gibbs_iter, mh_proposal_sd, bgp_max);
        if (i % thinning == 0) {
            cout << "Iter: " << i << endl;
            cout << "Curr params: " << params.get_bgp() << endl;
            cout << "Curr state: " << state->to_string() << endl;
            cout << "Curr log_lik: " << state->get_log_lik() << endl;
            
            if (i > burn_in) {
                bgps.push_back(params.get_bgp());
                fbps.push_back(params.get_fbp());

                // store state
                for (size_t g = 0; g < n_genes; g++) {
                    gsl_matrix_set(states, sample_idx, g, state->get_pathway_membership_of(g));
                }
                sample_idx++;
            }
        }
        if (i < burn_in) {
            bgps.push_back(params.get_bgp());
            fbps.push_back(params.get_fbp());
        }

        if (state->get_log_lik() > best_log_lik) {
            best_log_lik = state->get_log_lik();
            cout << "=====" << endl;
            cout << state->to_string() << endl;
            cout << state->get_log_lik() << endl;
            cout << "=====" << endl;
        }
        
        // perform simple adaptation during burn-in
        if (i < burn_in) {
            mh_proposal_sd = params.get_bgp()/2;
        } else if (i == burn_in) {
            // compute empirical standard deviation of the parameters sampled and use it
            mh_proposal_sd = gsl_stats_sd(bgps.data(), 1, bgps.size());
            bgps.clear();
            fbps.clear();
        }
    }
    
    // store the last state
    bgps.push_back(params.get_bgp());
    fbps.push_back(params.get_fbp());
    
    // store state
    for (size_t g = 0; g < n_genes; g++) {
        gsl_matrix_set(states, sample_idx, g, state->get_pathway_membership_of(g));
    }
    sample_idx++;

    
    assert(sample_idx == states->size1);

    delete [] pathway_indices;
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
            double mh_proposal_sd,
            bool use_lik_tempering,
            unsigned int n_threads)
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

    run_pg_from_matrix(seed, obs_matrix, model_len, n_pg_iter, n_particles, n_smc_iter, n_kernel_iter, n_mh_w_gibbs_iter, has_passenger, swap_prob, fbp_max, bgp_max, 0, 0, mh_proposal_sd, use_lik_tempering, output_path, n_threads);

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
                                                                                      bool use_lik_tempering,
                                                                                      const char *output_path,
                                                                                      unsigned int n_threads)
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
    smc_options.csmc_set_partile_population = true;
    smc_options.main_seed = gsl_rng_get(random);
    smc_options.resampling_seed = gsl_rng_get(random);
    smc_options.num_threads = n_threads;

    // pmcmc options
    PMCMCOptions pmcmc_options(gsl_rng_get(random), n_pg_iter);

    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, use_lik_tempering);

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
    gsl_rng_free(random);
    return pg;
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
                bool use_lik_tempering,
                unsigned int n_threads)
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

    return run_smc_from_matrix(seed, obs_matrix, model_len, n_particles, n_smc_iter, n_kernel_iter, has_passenger, swap_prob, fbp, bgp, states, log_weights, use_lik_tempering, n_threads);
    
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
                           bool use_lik_tempering,
                           unsigned int n_threads)
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
    smc_options.num_threads = n_threads;
    smc_options.debug = true;

    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, use_lik_tempering);

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

    gsl_rng_free(random);

    return smc.get_log_marginal_likelihood();
}

void process_samples(const string output_path,
                     const gsl_matrix &data,
                     const gsl_matrix *states,
                     vector<double> &posterior_bgps,
                     vector<double> &posterior_fbps,
                     unsigned int n_drivers,
                     bool has_passenger,
                     double prior_passenger_prob)
{
    unsigned int n_genes = data.size2;
    unsigned int n_states = states->size1;

    // find posterior mean
    vector<double> fbp_mean(1), bgp_mean(1);
    fbp_mean[0] = gsl_stats_mean(posterior_fbps.data(), 1, posterior_fbps.size());
    bgp_mean[0] = gsl_stats_mean(posterior_bgps.data(), 1, posterior_bgps.size());

    // compute log_lik for each of the states 1) at posterior mean and 2) theta_i
    // compute log(model_evidence)
    unsigned int *pathway = new unsigned int[n_genes];
    vector<double> log_liks_at_posterior_mean(n_states);
    double log_model_evidence = DOUBLE_NEG_INF;
    double log_lik_at_i;
    //double log_prior_at_i;
    for (size_t i = 0; i < n_states; i++) {
        for (size_t n = 0; n < n_genes; n++) {
            pathway[n] = (unsigned int)gsl_matrix_get(states, i, n);
        }
        log_liks_at_posterior_mean[i] = compute_likelihood_from_matrix(&data, pathway, n_drivers, n_genes, has_passenger, fbp_mean[0], bgp_mean[0]);
        log_lik_at_i = compute_likelihood_from_matrix(&data, pathway, n_drivers, n_genes, has_passenger, posterior_fbps[i], posterior_bgps[i]);
        log_model_evidence = log_add(log_model_evidence, -log_lik_at_i);
//        log_prior_at_i = log_prior(n_drivers, n_genes, has_passenger, prior_passenger_prob, pathway);
        //log_model_evidence = log_add(log_model_evidence, -(log_lik_at_i + log_prior_at_i));
    }
    log_model_evidence -= log(n_states);
    log_model_evidence *= -1;
    cout << "Model len: " << n_drivers << ". log p(Y): " << log_model_evidence << endl;
    
    // output posterior mean
    // output likelihood at posterior mean
    // output model evidence estimate
    write_vector(output_path + "/posterior_fbp.csv", fbp_mean);
    write_vector(output_path + "/posterior_bgp.csv", bgp_mean);
    write_vector(output_path + "/log_liks_at_posterior_mean.csv", log_liks_at_posterior_mean);
    vector<double> vec;
    vec.push_back(log_model_evidence);
    write_vector(output_path + "/log_model_evidence.csv", vec);
    
    delete [] pathway;
}

void perform_prediction(const char *mcmc_samples_path,
                        const char *data_path,
                        unsigned int model_len,
                        bool has_passenger,
                        unsigned int *true_pathway)
{
    // read the file containing the pathway samples
    string sample_path_str(mcmc_samples_path, strlen(mcmc_samples_path));
    vector<vector<unsigned int> > states = read_states(sample_path_str + "/states.csv", ",");
    cout << states.size() << ", " << states[0].size() << endl;

    // read data for the patients to be predicted
    string data_path_str(data_path, strlen(data_path));
    gsl_matrix *data = read_data(data_path, ",");
    size_t n_patients = data->size1;
    size_t n_genes = states[0].size();
    cout << "Prediction on " << n_patients << " patients." << endl;

    // read posterior samples for the parameters
    vector<double> fbps, bgps;
    read_error_params(sample_path_str + "/fbps.csv", fbps);
    read_error_params(sample_path_str + "/bgps.csv", bgps);

    // compute the row sum for each patient
    vector<unsigned int> row_sum(n_patients);
    compute_row_sum(*data, row_sum);

    vector<LinearProgressionState> lpm_states;
    for (size_t i = 0; i < states.size(); i++) {
        LinearProgressionState state(*data, row_sum, n_genes, model_len, has_passenger);
        for (unsigned int n = 0; n < n_genes; n++) {
            state.update_pathway_membership(n, states[i][n]);
        }
        lpm_states.push_back(state);
    }
    assert(lpm_states.size() == bgps.size());
    
    LinearProgressionState true_state(*data, row_sum, n_genes, model_len, has_passenger);
    if (true_pathway != 0) {
        for (unsigned int n = 0; n < n_genes; n++) {
            true_state.update_pathway_membership(n, true_pathway[n]);
        }
    }
    
    double log_sum = DOUBLE_NEG_INF;
    double log_val;
    for (size_t m = 0; m < n_patients; m++) {
        double max = DOUBLE_NEG_INF;
        size_t best_stage = 0;
        for (size_t stage = 0; stage < model_len; stage++) {
            log_sum = DOUBLE_NEG_INF;
            for (size_t i = 0; i < lpm_states.size(); i++) {
                log_val = compute_likelihood_for_sample(lpm_states[i].get_cache_at(m), lpm_states[i], stage, bgps[i], fbps[i]);
                log_sum = log_add(log_sum, log_val);
                
                if (true_pathway != 0) {
                    log_val = compute_likelihood_for_sample(true_state.get_cache_at(m), true_state, stage, bgps[i], fbps[i]);
                }
            }
            if (log_sum > max) {
                max = log_sum;
                best_stage = stage;
            }
        }
    }
    cout << endl;
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
                       unsigned int n_mc_jobs,
                       unsigned int n_smc_threads,
                       double *log_marginal_sum,
                       double *log_marginal_smc,
                       bool use_lik_tempering)
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
    
    // generate seeds for SMC
    vector<long> main_seeds, resampling_seeds;
    for (unsigned int i = 0; i < n_mc_samples; i++) {
        main_seeds.push_back(gsl_rng_get(random));
        resampling_seeds.push_back(gsl_rng_get(random));
    }
    
    //double log_prior = log_pathway_uniform_prior(model_len, n_genes);
    vector<double> log_marginals(n_mc_samples);
    
    SMCOptions smc_options;
    smc_options.num_particles = n_particles;
    smc_options.ess_threshold = 1.0;
    smc_options.resample_last_round = false;
    smc_options.debug = false;
    smc_options.csmc_set_partile_population = false;
    smc_options.main_seed = 1;
    smc_options.resampling_seed = 2;
    smc_options.num_threads = n_smc_threads;

    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, use_lik_tempering);

    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
    LinearProgressionParameters params(0.05, 0.05);
    shared_ptr<ParticleGenealogy<LinearProgressionState> > genealogy = csmc.initialize(params);
    cout << csmc.get_log_marginal_likelihood() << endl;
    for (size_t i = 0; i < 100; i++) {
        genealogy = csmc.run_csmc(params, genealogy);
        cout << csmc.get_log_marginal_likelihood() << endl;
        const LinearProgressionState * state = genealogy.get()->get_state_ptr_at(n_smc_iter - 1).get();
        cout << state->to_string() << endl;
    }

//    auto start = std::chrono::high_resolution_clock::now();
//
//    printf("n_unique_states, fbp, bgp, log_marginal_sum, log_marginal_smc\n");
//    omp_set_num_threads(n_mc_jobs);
//#pragma omp parallel for
//    for (unsigned int n = 0; n < n_mc_samples; n++) {
//        // smc options
//        SMCOptions smc_options;
//        smc_options.num_particles = n_particles;
//        smc_options.ess_threshold = 1.0;
//        smc_options.resample_last_round = false;
//        smc_options.debug = false;
//        smc_options.csmc_set_partile_population = false;
//        smc_options.main_seed = main_seeds[n];
//        smc_options.resampling_seed = resampling_seeds[n];
//        smc_options.num_threads = n_smc_threads;
//
//        // LPM model proposal
//        LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, use_lik_tempering);
//
//        // Declare conditional SMC object
//        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
//
//        // Construct LPMParam object
//        LinearProgressionParameters param(fbps[n], bgps[n]);
//
//        // run SMC
//        csmc.initialize(param);
//
//        // get log marginal estimate from SMC
//        if (log_marginal_smc != 0) {
//            log_marginal_smc[n] = csmc.get_log_marginal_likelihood();
//        }
//
//        // enumerate over all particles that was ever generated to approximate \sum_x p(y | x, \theta)
//        unordered_set<string> unique_states;
//        log_marginals[n] = DOUBLE_NEG_INF;
//        for (unsigned int r = 0; r < n_smc_iter; r++) {
//            for (unsigned int i = 0; i < n_particles; i++) {
//                const LinearProgressionState &state = csmc.get_state(r, i);
//                string key = state.to_string();
//                if (unique_states.count(key) == 0) {
//                    log_marginals[n] = log_add(log_marginals[n], state.get_log_lik());
//                    unique_states.insert(key);
//                }
//            }
//        }
//
//        log_marginals[n] += log_prior;
//        if (log_marginal_sum != 0) {
//            log_marginal_sum[n] = log_marginals[n];
//        }
//    }
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    cout << elapsed.count() << " seconds elapsed." <<  endl;

    // combine log_marginals to approximate:
    // log(fhat) = \int_{\theta} \sum_x p(y|x, k, \theta) p(x|k) p(\theta) d(\theta)
//    double log_f_hat = DOUBLE_NEG_INF;
//    for (unsigned int n = 0; n < n_mc_samples; n++) {
//        log_f_hat = log_add(log_f_hat, log_marginals[n]);
//    }
    
    gsl_rng_free(random);

    //return log_f_hat;
    return 0.0;
    
    delete obs_matrix;
}

double model_selection(long seed,
                       const char *dat_file,
                       const char *output_path,
                       unsigned int model_len,
                       unsigned int n_mc_samples,
                       unsigned int n_particles,
                       unsigned int n_smc_iter,
                       unsigned int n_kernel_iter,
                       bool has_passenger,
                       double swap_prob,
                       const double *fbps,
                       const double *bgps,
                       unsigned int n_mc_jobs,
                       unsigned int n_smc_threads,
                       bool use_lik_tempering)
{
    double *log_marginal_sum = new double[n_mc_samples];
    double *log_marginal_smc = new double[n_mc_samples];
    double fhat = model_selection(seed, dat_file, model_len, n_mc_samples, n_particles, n_smc_iter, n_kernel_iter, has_passenger, swap_prob, fbps, bgps, n_mc_jobs, n_smc_threads, log_marginal_sum, log_marginal_smc);
    
    // output log_marginal_sum, log_marginal_smc
    // model_len, fbps[i], bgps[i], lgo_marginal_sum[i], log_marginal_smc[i]
    gsl_matrix *mat = gsl_matrix_alloc(n_mc_samples, 5);
    for (size_t i = 0; i < n_mc_samples; i++) {
        gsl_matrix_set(mat, i, 0, model_len);
        gsl_matrix_set(mat, i, 1, fbps[i]);
        gsl_matrix_set(mat, i, 2, bgps[i]);
        gsl_matrix_set(mat, i, 3, log_marginal_sum[i]);
        gsl_matrix_set(mat, i, 4, log_marginal_smc[i]);
    }
    write_matrix_as_csv(output_path, *mat);
    gsl_matrix_free(mat);

    return fhat;
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

double compute_likelihood_from_matrix(const gsl_matrix *obs_matrix,
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
    string stages_file = output_path_str + "/patient_stages.csv";
    string pathway_file = output_path_str + "/generative_mem_mat.csv";
    string clean_matrix_file = output_path_str + "/clean_matrix.csv";
    string data_matrix_file = output_path_str + "/matrix.csv";

    write_vector(stages_file, stages);
    write_vector(pathway_file, true_pathway);
    write_matrix_as_csv(clean_matrix_file, *data_matrix);
    
    add_noise(random, fbp, bgp, data_matrix);

    // output: data after contamination
    write_matrix_as_csv(data_matrix_file, *data_matrix);
    
    gsl_rng_free(random);
    
    delete data_matrix;
}
