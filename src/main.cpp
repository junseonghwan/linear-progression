
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>

#include <spf/pmcmc.hpp>
#include <spf/particle_population.hpp>
#include <spf/sampling_utils.hpp>
#include <spf/csmc.hpp>
#include <spf/smc.hpp>
#include <spf/spf.hpp>

#include "lpm_likelihood.hpp"
#include "lpm_model.hpp"
#include "lpm_params.hpp"
#include "lpm_state.hpp"

using namespace std;

const double TOL = 1e-3;

double *read_data(string file_name, size_t n_genes, size_t n_patients);
vector<size_t> *read_stages(string file_name);
vector<size_t> *read_ground_truth(string path_name, size_t K, size_t N);

void exp_infer_pathway(gsl_rng *random)
{
    string file_name = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients320/rep1/dat.csv";
    string ground_truth_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients320/rep1/";
    string stages_file_name = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients320/rep1/stages.csv";
    size_t n_patients = 320;
    size_t n_genes = 25;
    size_t n_smc_iter = 100;
    size_t n_mh_iter = 1;
    size_t true_model_length = 5;
    size_t max_model_length = 8;
    double fbp = 0.1;
    double bgp = 0.1;

    double *obs = read_data(file_name, n_genes, n_patients);
    gsl_matrix_view obs_matrix = gsl_matrix_view_array(obs, n_patients, n_genes);
    vector<size_t> *stages = read_stages(stages_file_name);
    vector<size_t> *gt = read_ground_truth(ground_truth_path, true_model_length, n_genes);
    LinearProgressionState true_state(n_genes, true_model_length, false);
    true_state.initialize(random);
    for (size_t g = 0; g < gt->size(); g++) {
        true_state.update_pathway_membership(g, gt->at(g));
    }
    
    LinearProgressionParameters params(fbp, bgp, *stages);
    
    double log_lik_at_truth = compute_pathway_likelihood(obs_matrix, true_state, params);
    cout << "log lik at truth: " << log_lik_at_truth << endl;

    SMCOptions smc_options;
    smc_options.ess_threshold = 1.0;
    smc_options.num_particles = 200;
    smc_options.resample_last_round = false;

    LinearProgressionState *best_state = 0;
    for (size_t l = 1; l <= max_model_length; l++) {
        cout << "==== Model length: " << l << "====" << endl;
        LinearProgressionModel *model = new LinearProgressionModel(n_genes, l, n_smc_iter, n_mh_iter, &obs_matrix, false);
//        if (l == 1) {
//            model = new LinearProgressionModel(n_genes, l, n_smc_iter, n_mh_iter, &obs_matrix, false);
//        } else {
//            model = new LinearProgressionModel(n_genes, l, n_smc_iter, n_mh_iter, &obs_matrix, best_state, false);
//        }
        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(model, &smc_options);
        ParticleGenealogy<LinearProgressionState> *genealogy = csmc.initialize(params);
        vector<LinearProgressionState> *particles = csmc.get_particles(n_smc_iter - 1);
        size_t max_iter = l == 1 ? 1 : 4;
        vector<double> *log_weights = csmc.get_log_weights(n_smc_iter - 1);
        double logZ = csmc.get_log_marginal_likelihood();
        cout << logZ << endl;
        auto *pop = new ParticlePopulation<LinearProgressionState>(particles, log_weights);
        for (size_t num_iter = 0; num_iter < max_iter; num_iter++) {
            model->set_initial_population(pop);
            genealogy = csmc.run_csmc(params, genealogy);
            logZ = csmc.get_log_marginal_likelihood();
            cout << logZ << endl;
            particles = csmc.get_particles(n_smc_iter - 1);
            log_weights = csmc.get_log_weights(n_smc_iter - 1);
            pop = new ParticlePopulation<LinearProgressionState>(particles, log_weights);
        }

        particles = csmc.get_particles(n_smc_iter - 1);
        log_weights = csmc.get_log_weights(n_smc_iter - 1);
        double best_logw = DOUBLE_NEG_INF;
        // find the best state
        for (size_t i = 0; i < particles->size(); i++) {
            LinearProgressionState &state = particles->at(i);
            if (log_weights->at(i) > best_logw) {
                best_logw = log_weights->at(i);
                best_state = &state;
                cout << log_weights->at(i) << endl;
                cout << state.to_string() << endl;
            }
        }
    }
}

void exp_infer_model_length(gsl_rng *random)
{
    string file_name = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients320/rep0/dat.csv";
    //string stages_file_name = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/LPM/simulated_data/model_selection/npatients320/rep0/stages.csv";
    size_t n_patients = 320;
    size_t n_genes = 25;
    size_t n_smc_iter = 5;
    //size_t true_model_length = 5;
    
    double *obs = read_data(file_name, n_genes, n_patients);
    gsl_matrix_view obs_matrix = gsl_matrix_view_array(obs, n_patients, n_genes);
    //vector<size_t> *stages = read_stages(stages_file_name);

    SMCOptions smc_options;
    smc_options.ess_threshold = 1.0;
    smc_options.num_particles = 100;
    smc_options.resample_last_round = false;

    // sample \delta, \epsilon from prior say, uniform on (0, 0.2) x (0, 0.2)
    size_t num_samples = 20;
    //LinearProgressionParameters params(0.1, 0.1, *stages);
    LinearProgressionParameters params(0.1, 0.1);
    for (size_t l = 1; l < 9; l++) {
        LinearProgressionModel *model = new LinearProgressionModel(n_genes, l, n_smc_iter, 10, &obs_matrix, false);
        double model_evidence = 0.0;
        for (size_t i = 0; i < num_samples; i++) {
            double fbp = gsl_ran_flat(random, 0.0, 0.2);
            double bgp = gsl_ran_flat(random, 0.0, 0.2);
            params.set_bgp(bgp);
            params.set_fbp(fbp);
            SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
            smc.run_smc(params);
            cout << smc.get_log_marginal_likelihood() << endl;
            model_evidence += smc.get_log_marginal_likelihood();
        }
        cout << "Model " << l << ": " << (model_evidence / num_samples) << endl;
    }
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

vector<size_t> *read_ground_truth(string path_name, size_t K, size_t N)
{
    vector<size_t> *gt = new vector<size_t>(N);
    for (size_t k = 1; k <= K; k++) {
        string line;
        string file_name = path_name + "/pathway" + to_string(k) + ".csv";
        ifstream dat_file (file_name);
        if (!dat_file.is_open())
        {
            cerr << "Could not open the file: " << file_name << endl;
            exit(-1);
        }
        
        while ( getline (dat_file, line) )
        {
            unsigned int gene = stoi(line);
            (*gt)[gene] = k - 1;
        }
        dat_file.close();
    }

    return gt;
}

bool test_likelihood1()
{
    double fbp = 0.072;
    double bgp = 0.0321;
    
    int n_patients = 2;
    int n_pathways = 3;
    int n_genes = 3;
    
    double y[] = {
        1, 0, 1,
        1, 1, 1,
    };
    
    double x[] = {
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    };
    
    double r[n_patients*n_pathways];
    int pathway_sizes[] = {1, 1, 1};
    
    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
    gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
    gsl_matrix_view R = gsl_matrix_view_array(r, n_patients, n_pathways);
    
    /* Compute C = A B */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, &Y.matrix, &X.matrix,
                    0.0, &R.matrix);
    
    
    // test likelihood calculation
    double truth[] = {-3.546249, -0.2241706};
    
    int stages[] = {0, 2};
    for (int m = 0; m < n_patients; m++) {
        double log_lik = 0.0;
        for (int k = 0; k < n_pathways; k++) {
            double val = gsl_matrix_get(&R.matrix, m, k);
            if (k <= stages[m]) {
                log_lik += compute_log_lik(val, pathway_sizes[k], bgp, fbp);
            } else {
                log_lik += compute_log_lik_bg(val, pathway_sizes[k], bgp, fbp);
            }
        }
        if (abs(truth[m] - log_lik) > TOL) {
            cerr << "Likelihood test 1: error in likelihood computation for sample " << m << endl;
            cerr << "Expected: " << truth[m] << endl;
            cerr << "Obtained: " << log_lik << endl;
            return false;
        }
    }
    cout << "Likelihood test 1 passed!" << endl;
    return true;
}

bool test_likelihood2()
{
    // simple test:
    // read in the data
    // read in the true pathway
    // compute the likelihood
    
    double fbp = 0.027;
    double bgp = 0.1321;
    
    int n_patients = 3;
    int n_pathways = 3;
    int n_genes = 7;
    
    double y[] = {
        1, 0, 1, 0, 0, 0, 0,
        1, 1, 1, 0, 1, 1, 1,
        0, 0, 1, 1, 0, 0, 0
    };
    
    double x[] = {
        1, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 0, 1
    };
    
    double r[n_patients*n_pathways];
    int pathway_sizes[] = {2, 4, 1};
    
    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
    gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
    gsl_matrix_view R = gsl_matrix_view_array(r, n_patients, n_pathways);
    
    /* Compute C = A B */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, &Y.matrix, &X.matrix,
                    0.0, &R.matrix);
    
    
    // test likelihood calculation
    double truth[] = {-4.21657, -4.503214, -4.839183};
    
    int stages[] = {2, 2, 1};
    for (int m = 0; m < n_patients; m++) {
        double log_lik = 0.0;
        for (int k = 0; k < n_pathways; k++) {
            double val = gsl_matrix_get(&R.matrix, m, k);
            if (k <= stages[m]) {
                log_lik += compute_log_lik(val, pathway_sizes[k], bgp, fbp);
            } else {
                log_lik += compute_log_lik_bg(val, pathway_sizes[k], bgp, fbp);
            }
        }
        if (abs(truth[m] - log_lik) > TOL) {
            cerr << "Likelihood test 2: error in likelihood computation for sample " << m << endl;
            cerr << "Expected: " << truth[m] << endl;
            cerr << "Obtained: " << log_lik << endl;
            return false;
        }
    }
    cout << "Likelihood test 2 passed!" << endl;
    return true;
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

bool test_lpm_model_importance_sampling()
{
    // for small number of genes and pathway, we can compute the exact value of the marginal likelihood by
    // enumerating over all possibilities
    // then, we can check that the marginal likelihood estimate from SMC approaches this value
    
    double fbp = 0.027;
    double bgp = 0.1321;
    
    int n_patients = 3;
    int n_pathways = 2;
    int n_genes = 3;
    
    double y[] = {
        1, 1, 0,
        1, 0, 1,
        1, 1, 1
    };
    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
    vector<size_t> stages(n_patients, 0);
    stages[0] = 1;
    stages[1] = 0;
    stages[2] = 2;
    
    // enumerate over possible values of x using 3 for loops
    double exact_log_lik = DOUBLE_NEG_INF;
    vector<size_t> pathway_membership(n_genes, 0);
    for (size_t i = 0; i < n_pathways; i++) {
        pathway_membership.at(0) = i;
        for (size_t j = 0; j < n_pathways; j++) {
            pathway_membership.at(1) = j;
            for (size_t k = 0; k < n_pathways; k++) {
                pathway_membership.at(2) = k;
                
                // convert pathway membership to matrix x
                double *x = convert_to_array(pathway_membership, n_pathways);
                gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
                double log_val = compute_pathway_likelihood(Y, X, stages, fbp, bgp);
                exact_log_lik = log_add(exact_log_lik, log_val);
                cout << "===" << endl;
                cout << log_val << ", " << exact_log_lik << endl;
                for (size_t l = 0; l < pathway_membership.size(); l++)
                    cout << pathway_membership[l] << endl;
                cout << "===" << endl;
            }
        }
    }
    exact_log_lik += log(1./8); // multiply by the prior
    cout << "Exact marginal likelihood: " << exact_log_lik << endl;
    
    // run importance sampling with proposal from the prior, to get an estimate of the marginal log likelihood
    SMCOptions smc_options;
    smc_options.ess_threshold = 1.0;
    smc_options.num_particles = 100000;
    smc_options.resample_last_round = false;
    
    size_t n_smc_iter = 1;
    bool allocate_passenger_pathway = false;
    LinearProgressionModel *model = new LinearProgressionModel(n_genes, n_pathways, n_smc_iter, 10, &Y, allocate_passenger_pathway);
    LinearProgressionParameters params(fbp, bgp);
    params.set_patient_progression_stages(stages);
    
    SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
    smc.run_smc(params);
    ParticlePopulation<LinearProgressionState> *pop = smc.get_curr_population();
    double logZ = pop->get_log_norm();
    cout << "Estimated logZ: " << logZ << endl;
    if (abs(exp(exact_log_lik) - exp(logZ)) > TOL) {
        cerr << "Importance sampling estimate is incorrect." << endl;
        return false;
    }
    return true;
}

bool test_lpm_model_smc()
{
    double fbp = 0.027;
    double bgp = 0.1321;
    
    int n_patients = 3;
    int n_pathways = 2;
    int n_genes = 3;
    
    double y[] = {
        1, 1, 0,
        1, 0, 1,
        1, 1, 1
    };
    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
    vector<size_t> stages(n_patients, 0);
    stages[0] = 1;
    stages[1] = 0;
    stages[2] = 2;
    
    // enumerate over possible values of x using 3 for loops
    double exact_log_lik = DOUBLE_NEG_INF;
    vector<size_t> pathway_membership(n_genes, 0);
    for (size_t i = 0; i < n_pathways; i++) {
        pathway_membership.at(0) = i;
        for (size_t j = 0; j < n_pathways; j++) {
            pathway_membership.at(1) = j;
            for (size_t k = 0; k < n_pathways; k++) {
                pathway_membership.at(2) = k;
                
                // convert pathway membership to matrix x
                double *x = convert_to_array(pathway_membership, n_pathways);
                gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
                double log_val = compute_pathway_likelihood(Y, X, stages, fbp, bgp);
                exact_log_lik = log_add(exact_log_lik, log_val);
                cout << "===" << endl;
                cout << log_val << ", " << exact_log_lik << endl;
                for (size_t l = 0; l < pathway_membership.size(); l++)
                    cout << pathway_membership[l] << endl;
                cout << "===" << endl;
            }
        }
    }
    exact_log_lik += log(1./8);
    cout << "Exact marginal likelihood: " << exact_log_lik << endl;
    
    // run importance sampling with proposal from the prior, to get an estimate of the marginal log likelihood
    SMCOptions smc_options;
    smc_options.ess_threshold = 0.8;
    smc_options.num_particles = 2000;
    smc_options.resample_last_round = false;
    smc_options.track_population = true;
    size_t n_smc_iter = 200;
    
    bool allocate_passenger_pathway = false;
    LinearProgressionModel *model = new LinearProgressionModel(n_genes, n_pathways, n_smc_iter, 5, &Y, allocate_passenger_pathway);
    LinearProgressionParameters params(fbp, bgp);
    params.set_patient_progression_stages(stages);
    
    SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
    smc.run_smc(params);
    double logZ = smc.get_log_marginal_likelihood();
    cout << "=====" << endl;
    cout << "Estimated logZ: " << logZ << endl;
    if (abs(exp(exact_log_lik) - exp(logZ)) > TOL) {
        cerr << "SMC estimate is incorrect." << endl;
        return false;
    }
    return true;
}

void run_tests()
{
    if (!test_likelihood1()) {
        cerr << "Error: likelihood test 1" << endl;
        exit(-1);
    }
    if (!test_likelihood2()) {
        cerr << "Error: likelihood test 2" << endl;
        exit(-1);
    }
    if (!test_lpm_model_importance_sampling()) {
        cerr << "Error: importance sampling" << endl;
        exit(-1);
    }
    if (!test_lpm_model_smc()) {
        cerr << "Error: smc sampler" << endl;
        exit(-1);
    }
}

int main(int argc, char *argv[]) // TODO: take arguments from command line
{
    gsl_rng *random = generate_random_object(1);
    exp_infer_model_length(random);
    //exp_infer_pathway(random);
    return 0;
}
