
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>

#include <spf/particle_population.hpp>
#include <spf/pg.hpp>
#include <spf/sampling_utils.hpp>
#include <spf/csmc.hpp>
#include <spf/smc.hpp>
#include <spf/spf.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_model.hpp"
#include "lpm_params.hpp"
#include "lpm_pg_proposal.hpp"
#include "lpm_state.hpp"

using namespace std;

void compute_row_sum(gsl_matrix_view &obs, vector<size_t> &row_sum);
double *read_data(string file_name, size_t n_genes, size_t n_patients);
vector<size_t> *read_stages(string file_name);
vector<size_t> *read_ground_truth(string path_name);
void read_error_params(string path_name, double &fbp, double &bgp);

void exp_particle_gibbs(gsl_rng *random, string data_path, size_t n_patients, size_t n_genes, size_t true_model_length, size_t max_model_length)
{
    string obs_file = data_path + "/dat.csv";
    string ground_truth_file = data_path + "/pathways.csv";
    string stages_file = data_path + "/stages.csv";
    string params_file = data_path + "/params.csv";

    size_t n_smc_iter = 10;
    size_t n_mh_iter = 5;

    double *obs = read_data(obs_file, n_genes, n_patients);
    gsl_matrix_view obs_matrix_view = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> row_sum(n_patients);
    compute_row_sum(obs_matrix_view, row_sum);

    //vector<size_t> *true_stages = read_stages(stages_file);

    vector<size_t> *gt = read_ground_truth(ground_truth_file);
    LinearProgressionState true_state(n_genes, true_model_length, true);
    true_state.sample_from_prior(random);
    for (size_t g = 0; g < gt->size(); g++) {
        true_state.update_pathway_membership(g, gt->at(g));
    }
    double fbp = 0.0;
    double bgp = 0.0;
    read_error_params(params_file, fbp, bgp);

    // compute likelihood of the truth
    //LinearProgressionParameters true_params(fbp, bgp, *true_stages);
    LinearProgressionParameters true_params(fbp, bgp);
    double log_lik_at_truth = compute_pathway_likelihood(obs_matrix_view, row_sum, true_state, true_params);
    cout << log_lik_at_truth << endl;

    SMCOptions smc_options;
    smc_options.ess_threshold = 1.0;
    smc_options.num_particles = 100;
    smc_options.resample_last_round = false;

    PMCMCOptions pmcmc_options(1, 20);
    LinearProgressionState *best_state = 0;
    for (size_t l = 1; l <= max_model_length; l++) {
        LinearProgressionModel smc_model(n_genes, l, n_smc_iter, n_mh_iter, obs_matrix_view, row_sum, 0.2, true, false);
        if (l > 1) {
            smc_model.set_initial_state(best_state);
        }
        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(&smc_model, &smc_options);
        LPMParamProposal pg_proposal(obs_matrix_view, row_sum, n_patients, l, n_genes, 20);
        ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(&pmcmc_options, &csmc, &pg_proposal);
        pg.run();
        vector<ParticleGenealogy<LinearProgressionState> *> *states = pg.get_states();
        double best_log_lik = DOUBLE_NEG_INF;
        for (size_t i = 0; i < states->size(); i++) {
            ParticleGenealogy<LinearProgressionState> *genealogy = states->at(i);
            LinearProgressionState &state = genealogy->at(genealogy->size() - 1)->first;
            double log_lik = compute_pathway_likelihood(obs_matrix_view, row_sum, state, true_params);
            cout << log_lik << endl;
            cout << state.to_string() << endl;
            if (log_lik > best_log_lik) {
                best_log_lik = log_lik;
                best_state = &state;
            }
        }
    }
}

void run_particle_gibbs(gsl_rng *random, string data_path, size_t n_patients, size_t n_genes, size_t model_length)
{
    string obs_file = data_path + "/dat.csv";
    
    size_t n_smc_iter = 20;
    size_t n_mh_iter = 1;

    double *obs = read_data(obs_file, n_genes, n_patients);
    gsl_matrix_view obs_matrix = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> row_sum(n_patients);
    compute_row_sum(obs_matrix, row_sum);
    
    LinearProgressionModel smc_model(n_genes, model_length, n_smc_iter, n_mh_iter, obs_matrix, row_sum, 0.2, true, false);
    SMCOptions smc_options;
    smc_options.ess_threshold = 1.0;
    smc_options.num_particles = 200;
    smc_options.resample_last_round = false;
    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(&smc_model, &smc_options);

    PMCMCOptions pmcmc_options(1, 100);
    LPMParamProposal pg_proposal(obs_matrix, row_sum, n_patients, model_length, n_genes, 20);
    ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(&pmcmc_options, &csmc, &pg_proposal);
    pg.run();
    vector<ParticleGenealogy<LinearProgressionState> *> *states = pg.get_states();
    // marginalize over the states
    //vector<size_t> stages = pg_proposal.get_stages();
    for (size_t i = 0; i < states->size(); i++) {
        ParticleGenealogy<LinearProgressionState> *genealogy = states->at(i);
        LinearProgressionState &state = genealogy->at(genealogy->size() - 1)->first;
        cout << state.to_string() << endl;
    }
}

void exp_infer_pathway(gsl_rng *random, string data_path, size_t n_patients, size_t n_genes, size_t true_model_length)
{
    string file_name = data_path + "/dat.csv";
    string ground_truth_path = data_path + "/pathways.csv";
    string stages_file_name = data_path + "/stages.csv";

    size_t n_smc_iter = 10;
    size_t n_mh_iter = 1;
    size_t max_model_length = 8;
    double fbp = 0.17233289408058203;
    double bgp = 0.1685351284814174;

    double *obs = read_data(file_name, n_genes, n_patients);
    gsl_matrix_view obs_matrix = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> row_sum(n_patients);
    compute_row_sum(obs_matrix, row_sum);
    vector<size_t> *stages = read_stages(stages_file_name);
    vector<size_t> *gt = read_ground_truth(ground_truth_path);
    LinearProgressionState true_state(n_genes, true_model_length, true);
    true_state.sample_from_prior(random);
    for (size_t g = 0; g < gt->size(); g++) {
        true_state.update_pathway_membership(g, gt->at(g));
    }

    LinearProgressionParameters params(fbp, bgp, *stages);

    double log_lik_at_truth = compute_pathway_likelihood(obs_matrix, row_sum, true_state, params);
    cout << "log lik at truth: " << log_lik_at_truth << endl;

    SMCOptions smc_options;
    smc_options.ess_threshold = 1.0;
    smc_options.num_particles = 200;
    smc_options.resample_last_round = false;

    for (size_t l = 1; l <= max_model_length; l++) {
        cout << "==== Model length: " << l << "====" << endl;
        LinearProgressionModel *model = new LinearProgressionModel(n_genes, l, n_smc_iter, n_mh_iter, obs_matrix, row_sum, 0.2, true, true);
        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(model, &smc_options);
        ParticleGenealogy<LinearProgressionState> *genealogy = csmc.initialize(params);
        vector<LinearProgressionState> *particles = csmc.get_particles(n_smc_iter - 1);
        size_t max_iter = 5;
        vector<double> *log_weights = csmc.get_log_weights(n_smc_iter - 1);
        double logZ = csmc.get_log_marginal_likelihood();
        cout << logZ << endl;
        auto *pop = new ParticlePopulation<LinearProgressionState>(particles, log_weights);
        for (size_t num_iter = 0; num_iter < max_iter; num_iter++) {
            //model->set_initial_population(pop);
            genealogy = csmc.run_csmc(params, genealogy);
            logZ = csmc.get_log_marginal_likelihood();
            cout << logZ << endl;
            particles = csmc.get_particles(n_smc_iter - 1);
            log_weights = csmc.get_log_weights(n_smc_iter - 1);
            pop = new ParticlePopulation<LinearProgressionState>(particles, log_weights);
        }

        particles = csmc.get_particles(n_smc_iter - 1);
        log_weights = csmc.get_log_weights(n_smc_iter - 1);

        double best_log_lik = DOUBLE_NEG_INF;
        LinearProgressionState *best_state = 0;

        // find the best state
        for (size_t i = 0; i < particles->size(); i++) {
            LinearProgressionState &state = particles->at(i);
            double log_lik = state.get_log_lik();
            if (log_lik > best_log_lik) {
                best_log_lik = log_lik;
                best_state = &state;
                cout << log_lik << endl;
                cout << state.to_string() << endl;
            }
        }
    }
}

void exp_infer_model_length(gsl_rng *random, string data_path, size_t n_patients, size_t n_genes, size_t true_model_length)
{
    string file_name = data_path + "/dat.csv";
    string stages_file_name = data_path + "/stages.csv";
    string ground_truth_path = data_path + "/pathways.csv";

    size_t n_smc_iter = 10;
    size_t n_mcmc_iter = 5;

    double *obs = read_data(file_name, n_genes, n_patients);
    gsl_matrix_view obs_matrix = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> row_sum(n_patients);
    compute_row_sum(obs_matrix, row_sum);
    vector<size_t> *stages = read_stages(stages_file_name);

    vector<size_t> *gt = read_ground_truth(ground_truth_path);
    LinearProgressionState true_state(n_genes, true_model_length, false);
    true_state.sample_from_prior(random);
    for (size_t g = 0; g < gt->size(); g++) {
        true_state.update_pathway_membership(g, gt->at(g));
    }

    SMCOptions smc_options;
    smc_options.ess_threshold = 0.5;
    smc_options.num_particles = 50;
    smc_options.resample_last_round = false;

    // sample \delta, \epsilon from prior say, uniform on (0.05, 0.15) x (0.05, 0.15)
    size_t num_samples = 40;
    LinearProgressionParameters params(0.1, 0.1, *stages);
    vector<double> log_marginal(num_samples);
    cout << "model, fbp, bgp, logZ" << endl;
    for (size_t l = 1; l <= 9; l++) {
        LinearProgressionModel *model = new LinearProgressionModel(n_genes, l, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum);
        SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
        double model_evidence = DOUBLE_NEG_INF;
        for (size_t i = 0; i < num_samples; i++) {
            double fbp = gsl_ran_flat(random, 0, 0.3);
            double bgp = gsl_ran_flat(random, 0, 0.3);
            params.set_bgp(bgp);
            params.set_fbp(fbp);
            
            if (l == true_model_length) {
                double log_lik_at_truth = compute_pathway_likelihood(obs_matrix, row_sum, true_state, params);
                cout << true_state.to_string() << endl;
                cout << log_lik_at_truth << endl;
            }

            smc.run_smc(params);
            // find the best state and print its likelihood
            vector<LinearProgressionState> *particles = smc.get_curr_population()->get_particles();

            double best_logw = DOUBLE_NEG_INF;
            LinearProgressionState *best_state = 0;

            // find the best state
            for (size_t i = 0; i < particles->size(); i++) {
                LinearProgressionState &state = particles->at(i);
                if (state.get_log_lik() > best_logw) {
                    best_logw = state.get_log_lik();
                    best_state = &state;
                }
            }
            double best_state_log_lik = compute_pathway_likelihood(obs_matrix, row_sum, *best_state, params);
            cout << best_state->to_string() << endl;
            if (l == true_model_length) {
                cout << l << ", " << fbp << ", " << bgp << ", " << smc.get_log_marginal_likelihood() << ", " << best_logw << endl;
            } else {
                cout << l << ", " << fbp << ", " << bgp << ", " << smc.get_log_marginal_likelihood() << ", " << best_state_log_lik << endl;
            }
            model_evidence = log_add(model_evidence, smc.get_log_marginal_likelihood());
            log_marginal[i] = smc.get_log_marginal_likelihood();
        }
        double log_mean = model_evidence - log(num_samples);
        cout << "Model " << l << ": " << log_mean << endl;
        // compute the MC standard error?
//        double ssq = DOUBLE_NEG_INF;
//        for (size_t i = 0; i < num_samples; i++) {
//            double val = 2 * log_subtract(log_marginal[i], log_mean); // 2 * log(x_i - xmean)
//            ssq = log_add(ssq, val);
//        }
//        double log_var = ssq - log(num_samples - 1);
//        cout << "Model " << l << ": " << log_mean << ", " << log_var << endl;
    }
}

void infer_model_length(gsl_rng *random, string data_path, size_t n_patients, size_t n_genes, size_t max_L)
{
    if (max_L > n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }

    string file_name = data_path + "/dat.csv";
    
    size_t n_smc_iter = 10;
    size_t n_mcmc_iter = 5;
    
    double *obs = read_data(file_name, n_genes, n_patients);
    gsl_matrix_view obs_matrix = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> row_sum(n_patients);
    compute_row_sum(obs_matrix, row_sum);
    
    SMCOptions smc_options;
    smc_options.ess_threshold = 0.5;
    smc_options.num_particles = 50;
    smc_options.resample_last_round = false;
    
    // sample \delta, \epsilon from prior say, uniform on (0.05, 0.15) x (0.05, 0.15)
    size_t num_samples = 40;
    LinearProgressionParameters params(0.1, 0.1);
    vector<double> log_marginal(num_samples);
    cout << "model, fbp, bgp, logZ" << endl;
    for (size_t l = 1; l <= max_L; l++) {
        LinearProgressionModel *model = new LinearProgressionModel(n_genes, l, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum);
        SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
        double model_evidence = DOUBLE_NEG_INF;
        for (size_t i = 0; i < num_samples; i++) {
            double fbp = gsl_ran_flat(random, 0, 0.3);
            double bgp = gsl_ran_flat(random, 0, 0.3);
            params.set_bgp(bgp);
            params.set_fbp(fbp);
            smc.run_smc(params);
            // find the best state and print its likelihood
            vector<LinearProgressionState> *particles = smc.get_curr_population()->get_particles();
            
            double best_logw = DOUBLE_NEG_INF;
            LinearProgressionState *best_state = 0;
            
            // find the best state
            for (size_t i = 0; i < particles->size(); i++) {
                LinearProgressionState &state = particles->at(i);
                if (state.get_log_lik() > best_logw) {
                    best_logw = state.get_log_lik();
                    best_state = &state;
                }
            }
            double best_state_log_lik = compute_pathway_likelihood(obs_matrix, row_sum, *best_state, params);
            cout << best_state->to_string() << endl;
            cout << l << ", " << fbp << ", " << bgp << ", " << smc.get_log_marginal_likelihood() << ", " << best_state_log_lik << endl;
            model_evidence = log_add(model_evidence, smc.get_log_marginal_likelihood());
            log_marginal[i] = smc.get_log_marginal_likelihood();
        }
        double log_mean = model_evidence - log(num_samples);
        cout << "Model " << l << ": " << log_mean << endl;
        // compute the MC standard error?
        //        double ssq = DOUBLE_NEG_INF;
        //        for (size_t i = 0; i < num_samples; i++) {
        //            double val = 2 * log_subtract(log_marginal[i], log_mean); // 2 * log(x_i - xmean)
        //            ssq = log_add(ssq, val);
        //        }
        //        double log_var = ssq - log(num_samples - 1);
        //        cout << "Model " << l << ": " << log_mean << ", " << log_var << endl;
    }
}

void read_error_params(string path_name, double &fbp, double &bgp)
{
    string line;
    ifstream dat_file (path_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << path_name << endl;
    }
    getline (dat_file, line); // first line is just the header
    getline (dat_file, line);
    vector<string> results;
    boost::split(results, line, boost::is_any_of(","));
    fbp = stod(results[0]);
    bgp = stod(results[1]);
    dat_file.close();
}

// returns double array of length n_patients*n_gene
double *read_data(string file_name, size_t n_genes, size_t n_patients)
{
    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        return 0;
    }
    
    double *obs = new double[n_patients*n_genes];
    vector<string> results;
    unsigned int idx = 0;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(","));
        if (results.size() != n_genes) {
            cerr << "Error in the input file: " << file_name << endl;
            return 0;
        }
        for (size_t i = 0; i < n_genes; i++) {
            obs[idx] = stod(results[i]);
            idx++;
        }
    }
    dat_file.close();
    
    return obs;
}

vector<size_t> *read_stages(string file_name)
{
    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        return 0;
    }
    
    auto *stages = new vector<size_t>();
    while ( getline (dat_file, line) )
    {
        unsigned int stage = stoi(line);
        stages->push_back(stage);
    }
    dat_file.close();
    
    return stages;
}

vector<size_t> *read_ground_truth(string file_name)
{
    vector<size_t> *gt = new vector<size_t>();
    ifstream dat_file (file_name);
    string line;
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }
    vector<string> results;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(","));
        //unsigned int gene = stoi(results[0]);
        unsigned int pathway = stoi(results[1]) - 1;
        gt->push_back(pathway);
    }
    dat_file.close();

    return gt;
}

bool test_likelihood4()
{
    string data_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients20/rep0/";
    string obs_file_path = data_path + "/dat.csv";
    string stages_file_path = data_path + "/stages.csv";
    string gt_file_path = data_path + "/pathways.csv";
    string param_file_path = data_path + "/params.csv";
    
    size_t n_patients = 20;
    size_t n_genes = 25;
    size_t true_model_length = 5;
    
    double *obs = read_data(obs_file_path, n_genes, n_patients);
    gsl_matrix_view obs_matrix_view = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> row_sum(n_patients);
    compute_row_sum(obs_matrix_view, row_sum);
    vector<size_t> *stages = read_stages(stages_file_path);
    vector<size_t> *gt = read_ground_truth(gt_file_path);
    LinearProgressionState true_state(n_genes, true_model_length, false);
    for (size_t g = 0; g < gt->size(); g++) {
        true_state.update_pathway_membership(g, gt->at(g));
    }

    double fbp = 0.0, bgp = 0.0;
    read_error_params(param_file_path, fbp, bgp);
    LinearProgressionParameters true_params(fbp, bgp, *stages);
    // compute the log-likelihood of the data given the true stages
    double log_lik = compute_pathway_likelihood(obs_matrix_view, row_sum, true_state, true_params);
    double log_lik_expected = -256.2;
    double ret = abs(log_lik - log_lik_expected);
    cout << log_lik << " - " << log_lik_expected << "=" << ret << endl;
    if (ret > 1e-3) {
        return false;
    }
    
    LinearProgressionParameters true_params_no_stages(fbp, bgp);
    // compute the log-likelihood of the data by marginalizing out the stages
    double log_lik_marginal_stages = compute_pathway_likelihood(obs_matrix_view, row_sum, true_state, true_params_no_stages);
    double log_lik_marginal_stages_expected = -260.847;
    ret = abs(log_lik_marginal_stages - log_lik_marginal_stages_expected);
    cout << log_lik_marginal_stages << " - " << log_lik_marginal_stages_expected << " = " << ret << endl;
    if (ret > 1e-3) {
        return false;
    }
    return true;
}

void compute_row_sum(gsl_matrix_view &obs, vector<size_t> &row_sum)
{
    // compute the row sum for obs matrix
    size_t M = obs.matrix.size1;
    size_t N = obs.matrix.size2;
    if (row_sum.size() != M) {
        cerr << "Error: dimension for row_sum: " << row_sum.size() << " does not num rows of obs: " << M << endl;
        exit(-1);
    }
    for (size_t m = 0; m < M; m++) {
        size_t sum_m = 0;
        for (size_t n = 0; n < N; n++) {
            sum_m += (size_t)gsl_matrix_get(&obs.matrix, m, n);
        }
        row_sum[m] = sum_m;
    }
}

double *convert_to_array(vector<size_t> &pathway_membership, size_t n_pathways)
{
    size_t n_genes = pathway_membership.size();
    double *x = new double[n_genes * n_pathways];
    for (size_t n = 0; n < n_genes; n++) {
        size_t idx = pathway_membership[n];
        for (size_t k = 0; k < n_pathways; k++) {
            x[n * n_pathways + k] = idx == k ? 1 : 0;
        }
    }
    return x;
}

int main(int argc, char *argv[]) // TODO: take arguments from command line
{
    gsl_rng *random = generate_random_object(312312);
    //test_likelihood4();

    //size_t n_patients = 365;
    //size_t n_genes = 8;
    //size_t true_model_length = 5;
    //string data_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/data/";
    //infer_model_length(random, data_path, n_patients, n_genes, n_genes);
    //size_t inferred_model_length = 4;

    //exp_infer_pathway(random, data_path, n_patients, n_genes, true_model_length);

    size_t n_patients = 500;
    size_t n_genes = 25;
    size_t true_n_drivers = 4;
    string data_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients" + to_string(n_patients) + "/rep1/";
    //run_particle_gibbs(random, data_path, n_patients, n_genes, true_n_drivers);
    exp_particle_gibbs(random, data_path, n_patients, n_genes, true_n_drivers, true_n_drivers);
    //exp_infer_pathway(random, data_path, n_patients, n_genes, true_n_drivers);
    return 0;
}
