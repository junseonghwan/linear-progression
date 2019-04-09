//
//  lpm.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-28.
//

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>

#include <gsl/gsl_sf.h>

#include <omp.h>

#include <spf/pg.hpp>
#include <spf/sampling_utils.hpp>
#include <spf/smc.hpp>
#include <spf/csmc.hpp>

#include "data_util.hpp"
#include "lpm.hpp"
#include "lpm_likelihood.hpp"
#include "lpm_model.hpp"
#include "lpm_params.hpp"
#include "lpm_pg_proposal.hpp"
#include "lpm_state.hpp"

using namespace std;

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
            double bgp_max)
{
    // convert char* to string
    string dat_file_str(dat_file, strlen(dat_file));
    string output_path_str(output_path, strlen(output_path));
    string move_type_str = "MH"; // later, provide GIBBS as an option
    
    // load the data
    gsl_matrix *obs_matrix = read_data(dat_file_str, false);
    size_t n_patients = obs_matrix->size1;
    size_t n_genes = obs_matrix->size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }
    
    // output the settings
    cout << "Running PG {" << endl;
    cout << "\tInput data: " << dat_file_str << endl;
    cout << "\tOutput path: " << output_path_str << endl;
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
    vector<size_t> row_sum(n_patients);
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

    // set Markov kernel move type
    MoveType kernel_move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;
    
    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, kernel_move_type);
    
    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);

    // Declare param proposal object
    LPMParamProposal pg_proposal(*obs_matrix, row_sum, n_mh_w_gibbs_iter, fbp_max, bgp_max);
    
    // Initialize PG object
    ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(pmcmc_options, csmc, pg_proposal);
    
    // run (output timing in seconds)
    auto start = std::chrono::high_resolution_clock::now();
    pg.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds elapsed." <<  endl;

    // output the solution at the last iteration of PG
    vector<double> &log_marginal_likelihoods = pg.get_log_marginal_likelihoods();
    vector<shared_ptr<LinearProgressionParameters> > &params = pg.get_parameters();
    vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states = pg.get_states();
    ParticleGenealogy<LinearProgressionState> *genealogy = states.at(states.size() - 1).get();
    const LinearProgressionState &state = genealogy->get_state_at(genealogy->size() - 1);
    cout << "Current solution:\n" << state.to_string() << endl;

    // output the states and parameters
    write_pg_output(output_path_str, states, log_marginal_likelihoods, params, pg_proposal);
}

void run_pg_from_matrix(long seed,
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
                        vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &ret_states,
                        vector<shared_ptr<LinearProgressionParameters>> &ret_params)
{
    string move_type_str = "MH";
    size_t n_patients = obs_matrix->size1;
    size_t n_genes = obs_matrix->size2;
    
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
    vector<size_t> row_sum(n_patients);
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
    
    // set Markov kernel move type
    MoveType kernel_move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;
    
    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, kernel_move_type);
    
    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
    
    // Declare param proposal object
    LPMParamProposal pg_proposal(*obs_matrix, row_sum, n_mh_w_gibbs_iter, fbp_max, bgp_max);
    
    // Initialize PG object
    ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(pmcmc_options, csmc, pg_proposal);
    
    // run (output timing in seconds)
    auto start = std::chrono::high_resolution_clock::now();
    pg.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds elapsed." <<  endl;
    
    //vector<double> &log_marginal_likelihoods = pg.get_log_marginal_likelihoods();
    ret_params = pg.get_parameters();
    ret_states = pg.get_states();
    ParticleGenealogy<LinearProgressionState> *genealogy = ret_states.at(ret_states.size() - 1).get();
    const LinearProgressionState &state = genealogy->get_state_at(genealogy->size() - 1);
    cout << "Current solution:\n" << state.to_string() << endl;
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
                double *log_weights)
{
    // convert char* to string
    string dat_file_str(dat_file, strlen(dat_file));
    
    // load the data
    gsl_matrix *obs_matrix = read_data(dat_file_str, false);
    size_t n_patients = obs_matrix->size1;

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

    return run_smc_from_matrix(seed, obs_matrix, model_len, n_particles, n_smc_iter, n_kernel_iter, has_passenger, swap_prob, fbp, bgp, states, log_weights);
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
                           double *log_weights)
{
    size_t n_patients = obs_matrix->size1;
    size_t n_genes = obs_matrix->size2;

    if (model_len >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    // compute the row sum for each patient
    vector<size_t> row_sum(n_patients);
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

    // set Markov kernel move type
    string move_type_str = "MH"; // later, provide GIBBS as an option
    MoveType kernel_move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;

    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, kernel_move_type);

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

//    for (size_t r = 0; r < n_smc_iter; r++) {
//        cout << "log_norm, " << smc.get_log_norm(r) << endl;
//    }
//    for (size_t r = 0; r < n_smc_iter; r++) {
//        cout << "iter: " << r << endl;
//        unordered_map<string, unsigned int> unique_states;
//        unordered_map<string, double> unique_states_normalized_weights;
//        for (size_t k = 0; k < n_particles; k++) {
//            const LinearProgressionState &state = smc.get_state(r, k);
//            double prob = exp(smc.get_log_weight(r, k) - smc.get_log_norm(r));
//            if (unique_states.count(state.to_string()) == 0) {
//                unique_states_normalized_weights[state.to_string()] = prob;
//                unique_states[state.to_string()] = 0; // initialize the value to 0
//            }
//            unique_states[state.to_string()] += 1;
//            cout << state.to_string() << ", " << smc.get_log_weight(r, k) << ", " << prob <<  endl;
//        }
//        cout << unique_states.size() << endl;
//        double sum = 0.0;
//        for (auto it = unique_states.begin(); it != unique_states.end(); ++it) {
//            double prob = unique_states_normalized_weights[it->first] * it->second;
//            sum += prob;
//            cout << it->first << ", " << unique_states_normalized_weights[it->first] << ", " << it->second << endl;
//        }
//        cout << "sum: " << sum << endl;
//    }

    if (states != 0 && log_weights != 0) {
        // fill the states and log_weights
        // dimension of states should be n_particles x n_genes
        // dimension of the log_weights should be n_particles
        size_t idx = 0;
        for (size_t i = 0; i < n_particles; i++) {
            log_weights[i] = smc.get_log_weight(n_smc_iter - 1, i);
            const LinearProgressionState &state = smc.get_state(n_smc_iter - 1, i);
            const vector<size_t> pathway = state.get_pathway_membership();
            for (size_t j = 0; j < n_genes; j++) {
                idx = i*n_genes + j;
                states[idx] = pathway[j];
            }
        }
    }

    return smc.get_log_marginal_likelihood();
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
    string move_type_str = "MH"; // later, provide GIBBS as an option
    
    // load the data
    gsl_matrix *obs_matrix = read_data(dat_file_str, false);
    size_t n_patients = obs_matrix->size1;
    size_t n_genes = obs_matrix->size2;
    
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
    vector<size_t> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);
    
    // set Markov kernel move type
    MoveType kernel_move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;

    // allocate random object
    gsl_rng *random = generate_random_object(seed);
    
    // generate seeds and parameters
    vector<long> main_seeds, resampling_seeds;
    for (size_t i = 0; i < n_mc_samples; i++) {
        main_seeds.push_back(gsl_rng_get(random));
        resampling_seeds.push_back(gsl_rng_get(random));
    }
    
    vector<double> log_marginals(n_mc_samples);
    vector<double> smc_log_marginals(n_mc_samples);
    auto start = std::chrono::high_resolution_clock::now();
    
    printf("n_unique_states, fbp, bgp, log_marginal_sum, log_marginal_smc\n");
    omp_set_num_threads(n_threads);
#pragma omp parallel for
    for (size_t n = 0; n < n_mc_samples; n++) {
        // smc options
        SMCOptions smc_options;
        smc_options.num_particles = n_particles;
        smc_options.ess_threshold = 1.0;
        smc_options.resample_last_round = false;
        smc_options.main_seed = main_seeds[n];
        smc_options.resampling_seed = resampling_seeds[n];

        // LPM model proposal
        LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, kernel_move_type);

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
        for (size_t r = 0; r < n_smc_iter; r++) {
            for (size_t i = 0; i < n_particles; i++) {
                const LinearProgressionState &state = csmc.get_state(r, i);
                string key = state.to_string();
                if (unique_states.count(key) == 0) {
                    log_marginals[n] = log_add(log_marginals[n], state.get_log_lik());
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
    for (size_t n = 0; n < n_mc_samples; n++) {
        log_f_hat = log_add(log_f_hat, log_marginals[n]);
    }

    // account for p(x | model_len)
    log_f_hat += log_pathway_uniform_prior(model_len, n_genes);

    return log_f_hat;
    
    delete obs_matrix;
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
    gsl_matrix *obs_matrix = read_data(dat_file_str, false);
    return compute_likelihood_from_matrix(obs_matrix, pathway, model_len, n_genes, has_passenger, fbp, bgp);
}

double compute_likelihood_from_matrix(gsl_matrix *obs_matrix,
                                      unsigned int *pathway,
                                      unsigned int model_len,
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
 
    size_t n_patients = obs_matrix->size1;

    // compute the row sum for each patient
    vector<size_t> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);

    // construct LPState from the pathway
    LinearProgressionState state(*obs_matrix, row_sum, n_genes, model_len, has_passenger);
    for (size_t i = 0; i < n_genes; i++) {
        state.update_pathway_membership(i, pathway[i]);
    }

    // construct params
    LinearProgressionParameters params(fbp, bgp);

    double log_lik = compute_pathway_likelihood(*obs_matrix, row_sum, state, params);
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
    sample_pathway_uniform(random, model_len, true_pathway);
    gsl_matrix *data_matrix = simulate_data(random, model_len, true_pathway, stages);

    // output: stages, pathway, parameters, data matrix before contamination
    string stages_file = output_path_str + "/stages.csv";
    string pathway_file = output_path_str + "/true_pathways.csv";
    string clean_matrix_file = output_path_str + "/clean_matrix.csv";
    string data_matrix_file = output_path_str + "/data_matrix.csv";

    write_csv(stages_file, stages);
    write_csv(pathway_file, true_pathway);
    write_matrix(clean_matrix_file, *data_matrix);
    
    add_noise(random, fbp, bgp, data_matrix);

    // output: data after contamination
    write_matrix(data_matrix_file, *data_matrix);
}
