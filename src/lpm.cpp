//
//  lpm.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-28.
//

#include <stdio.h>
#include <string>

#include <spf/pg.hpp>
#include <spf/sampling_utils.hpp>
#include <spf/csmc.hpp>

#include "data_util.hpp"
#include "lpm.hpp"
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
    
    // output the settings
    cout << "Running PG {" << endl;
    cout << "\tInput data: " << dat_file_str << endl;
    cout << "\tOutput path: " << output_path_str << endl;
    cout << "\tModel length: " << model_len << endl;
    cout << "\tNum PG iter: " << n_pg_iter << endl;
    cout << "\tNum particles iter: " << n_particles << endl;
    cout << "\tNum SMC iter: " << n_smc_iter << endl;
    cout << "\tNum kernel iter: " << n_kernel_iter << endl;
    cout << "\tNum MHwGibbs iter: " << n_mh_w_gibbs_iter << endl;
    cout << "\tAllocate passenger pathway: " << to_string(has_passenger) << endl;
    cout << "\tSwap prob: " << swap_prob << endl;
    cout << "\tFBP max: " << fbp_max << endl;
    cout << "\tBGP max: " << bgp_max << endl;
    cout << "}" << endl;
    
    // load the data
    gsl_matrix *obs_matrix = read_data(dat_file_str, true);
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

    // pmcmc options
    PMCMCOptions pmcmc_options(gsl_rng_get(random), n_pg_iter);

    // set Markov kernel move type
    MoveType kernel_move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;
    
    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, kernel_move_type);
    
    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);

    // Declare param proposal object
    LPMParamProposal pg_proposal(n_patients, n_mh_w_gibbs_iter, fbp_max, bgp_max);
    
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
    string move_type_str = "MH"; // later, provide GIBBS as an option
    
    // output the settings
    cout << "Running SMC {" << endl;
    cout << "\tInput data: " << dat_file_str << endl;
    cout << "\tModel length: " << model_len << endl;
    cout << "\tNum particles iter: " << n_particles << endl;
    cout << "\tNum SMC iter: " << n_smc_iter << endl;
    cout << "\tNum kernel iter: " << n_kernel_iter << endl;
    cout << "\tAllocate passenger pathway: " << to_string(has_passenger) << endl;
    cout << "\tSwap prob: " << swap_prob << endl;
    cout << "\tFBP: " << fbp << endl;
    cout << "\tBGP: " << bgp << endl;
    cout << "}" << endl;
    
    // load the data
    gsl_matrix *obs_matrix = read_data(dat_file_str, true);
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
    MoveType kernel_move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;
    
    // LPM model proposal
    LinearProgressionModel smc_model(n_genes, model_len, n_smc_iter, n_kernel_iter, *obs_matrix, row_sum, swap_prob, has_passenger, kernel_move_type);
    
    // Declare conditional SMC object
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
    
    // Construct LPMParam object
    LinearProgressionParameters param(fbp, bgp);
    
    // run (output timing in seconds)
    auto start = std::chrono::high_resolution_clock::now();
    csmc.initialize(param);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds elapsed." <<  endl;
    
    // fill the states and log_weights
    // dimension of states should be n_particles x n_genes
    // dimension of the log_weights should be n_particles
    size_t idx = 0;
    for (size_t i = 0; i < n_particles; i++) {
        log_weights[i] = csmc.get_log_weight(n_smc_iter - 1, i);
        const LinearProgressionState &state = csmc.get_state(n_smc_iter - 1, i);
        const vector<size_t> pathway = state.get_pathway_membership();
        for (size_t j = 0; j < n_genes; j++) {
            idx = i*n_genes + j;
            states[idx] = pathway[j];
        }
    }
    
    return csmc.get_log_marginal_likelihood();
    
}
