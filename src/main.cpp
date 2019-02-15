#include <chrono>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <thread>
#include <unordered_map>
#include <vector>
#include <string>

#include <omp.h>

#include <boost/program_options.hpp>

#include <gsl/gsl_statistics.h>

#include <spf/particle_population.hpp>
#include <spf/pg.hpp>
#include <spf/sampling_utils.hpp>
#include <spf/csmc.hpp>
#include <spf/smc.hpp>
#include <spf/spf.hpp>

#include "data_util.hpp"
#include "lpm_likelihood.hpp"
#include "lpm_model.hpp"
#include "lpm_model_cache.hpp"
#include "lpm_params.hpp"
#include "lpm_pg_proposal.hpp"
#include "lpm_state.hpp"

using namespace std;

void compute_row_sum(gsl_matrix &obs, vector<size_t> &row_sum)
{
    // compute the row sum for obs matrix
    size_t M = obs.size1;
    size_t N = obs.size2;
    if (row_sum.size() != M) {
        cerr << "Error: dimension for row_sum: " << row_sum.size() << " does not num rows of obs: " << M << endl;
        exit(-1);
    }
    for (size_t m = 0; m < M; m++) {
        size_t sum_m = 0;
        for (size_t n = 0; n < N; n++) {
            sum_m += (size_t)gsl_matrix_get(&obs, m, n);
        }
        row_sum[m] = sum_m;
    }
}

void run_particle_gibbs(gsl_rng *random, gsl_matrix &obs_matrix, vector<size_t> &row_sum,
                        size_t target_model_length,
                        size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles,
                        size_t n_pmcmc_iter, size_t n_mh_gibbs_iter,
                        double fbp_max, double bgp_max,
                        double swap_prob, MoveType move_type,
                        string output_dir, size_t n_threads)
{
    size_t n_patients = obs_matrix.size1;
    size_t n_genes = obs_matrix.size2;

    if (target_model_length >= n_genes)
    {
        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
        exit(-1);
    }
    
    // smc settings
    SMCOptions smc_options;
    smc_options.num_particles = n_particles;
    smc_options.ess_threshold = 1.0;
    smc_options.resample_last_round = false;

    // pmcmc options
    size_t random_seed = gsl_rng_get(random);

    shared_ptr<LinearProgressionState> prev_state_ptr;
    for (size_t l = 1; l <= target_model_length; l++) {
        cout << "Running with model length " << l << endl;
        LinearProgressionModel smc_model(n_genes, l, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum, swap_prob, true, move_type);
        if (l > 1) {
            smc_model.set_initial_state(random, prev_state_ptr.get());
        }
        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
        // for l < target, we don't need to generate too many samples as convergence is not the concern
        PMCMCOptions pmcmc_options(random_seed, round(n_pmcmc_iter/10));
        LPMParamProposal pg_proposal(n_patients, n_mh_gibbs_iter, fbp_max, bgp_max);
        ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(pmcmc_options, csmc, pg_proposal);
        pg.run();

        // get the last genealogy, last state to initialize the next step
        vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states = pg.get_states();
        ParticleGenealogy<LinearProgressionState> *genealogy = states.at(states.size() - 1).get();
        const LinearProgressionState &state = genealogy->get_state_at(genealogy->size() - 1);
        prev_state_ptr.reset(new LinearProgressionState(state)); // make a copy
        cout << "Current solution:\n" << state.to_string() << endl;

        if (l == target_model_length) {
            vector<shared_ptr<LinearProgressionParameters> > &params = pg.get_parameters();
            write_pg_output(output_dir, states, params);
        }
    }
    // output parameters and states for the targe model length
    
}

size_t infer_model_length(gsl_rng *random, gsl_matrix &obs_matrix, vector<size_t> &row_sum,
                            size_t max_model_len, size_t num_samples,
                            size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles, size_t ess,
                            double fbp_max, double bgp_max,
                            double swap_prob, MoveType move_type,
                            string output_dir, size_t num_threads)
{
    size_t n_patients = obs_matrix.size1;
    size_t n_genes = obs_matrix.size2;
    
    cout << "Inferring model length: " << endl;
    cout << "# Patients: " << n_patients << endl;
    cout << "# Genes: " << n_genes << endl;
    cout << "# Particles: " << n_particles << endl;
    cout << "# MC Samples: " << num_samples << endl;
    cout << "# SMC iters: " << n_smc_iter << endl;
    cout << "# MCMC iters: " << n_mcmc_iter << endl;

    if (max_model_len >= n_genes)
    {
        cerr << "Maximum model length must be smaller than the number of genes." << endl;
        exit(-1);
    }

    SMCOptions smc_options;
    smc_options.ess_threshold = ess;
    smc_options.num_particles = n_particles;
    smc_options.resample_last_round = false;
    smc_options.debug = false;

    double *log_marginals = new double[num_samples];
    double best_marginal_log_lik = DOUBLE_NEG_INF;
    size_t best_l = 0;
    vector<double> log_marginal_mu(max_model_len);
    vector<double> log_marginal_sd(max_model_len);
    for (size_t l = 1; l <= max_model_len; l++) {
        cout << "Model length: " << l << endl;
        if (num_threads == 1) {
            cout << "Running times: {" << endl;
        }

        double model_evidence = DOUBLE_NEG_INF;

        omp_set_num_threads(num_threads);
#pragma omp parallel
        {
            size_t id = omp_get_thread_num();
            size_t n_threads = omp_get_num_threads();
            size_t start_idx = id * num_samples/n_threads;
            size_t end_idx = (id + 1) * num_samples/n_threads;
            if (id == n_threads - 1) {
                end_idx = num_samples;
            }
            //printf("%ld : [%ld, %ld)\n", id, start_idx, end_idx);
#pragma omp for
            for (size_t i = start_idx; i < end_idx; i++) {
                LinearProgressionModel model(n_genes, l, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum, swap_prob, true, move_type);
                LinearProgressionParameters params(gsl_ran_flat(random, 0.0, fbp_max), gsl_ran_flat(random, 0.0, bgp_max));
                ConditionalSMC<LinearProgressionState, LinearProgressionParameters> smc(model, smc_options);
                auto start = std::chrono::high_resolution_clock::now();
                smc.initialize(params);
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                if (num_threads == 1) {
                    printf("%f\n", elapsed.count());
                }
                
                log_marginals[i] = smc.get_log_marginal_likelihood();
#pragma omp critical
                model_evidence = log_add(model_evidence, log_marginals[i]);
            }
        }
        if (num_threads == 1) {
            cout << "}" << endl;
        }

        model_evidence = model_evidence - log(num_samples);
        cout << "Model evidence: " << model_evidence << endl;
        if (model_evidence > best_marginal_log_lik) {
            best_marginal_log_lik = model_evidence;
            best_l = l;
        }

        // compute mean and standard deviation log p(y|\theta^j)
        log_marginal_mu[l-1] = gsl_stats_mean(log_marginals, 1, num_samples);
        log_marginal_sd[l-1] = gsl_stats_sd_m(log_marginals, 1, num_samples, log_marginal_mu[l-1]);
    }
    
    string output_file = output_dir + "/model_selection.csv";
    write_model_selection_output_to_file(output_file, log_marginal_mu, log_marginal_sd);

    return best_l;
}

void run_exp(gsl_rng *random, string sim_dat_path,
             size_t max_model_len, size_t num_samples,
             size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles, size_t ess,
             size_t n_pmcmc_iter, size_t n_mh_w_gibbs_iter,
             double fbp_max, double bgp_max,
             double swap_prob, MoveType move_type,
             size_t n_threads)
{
    // path to files: data
    // ground truth: pathway, patient stages, params, model length used for generating the data
    string obs_file = sim_dat_path + "/dat.csv";
    string ground_truth_file = sim_dat_path + "/pathways.csv";
    string stages_file = sim_dat_path + "/stages.csv";
    string params_file = sim_dat_path + "/params.csv";
    string true_model_len_file = sim_dat_path + "/true_model_len.csv";

    // read the data
    gsl_matrix *obs_matrix = read_data(obs_file);
    size_t n_patients = obs_matrix->size1;
    size_t n_genes = obs_matrix->size2;

    // compute the row sum
    vector<size_t> row_sum(n_patients);
    compute_row_sum(*obs_matrix, row_sum);
    
    // read ground truth model length
    size_t true_model_length = read_true_model_length(true_model_len_file);
    
    // read the ground truth pathway
    vector<size_t> *gt = read_ground_truth(ground_truth_file);
    LinearProgressionState true_state(*obs_matrix, row_sum, n_genes, true_model_length, true);
    for (size_t g = 0; g < gt->size(); g++) {
        true_state.update_pathway_membership(g, gt->at(g));
    }
    
    // read the parameters used for generating the data
    double fbp = 0.0;
    double bgp = 0.0;
    read_error_params(params_file, fbp, bgp);

    // compute likelihood of the truth
    LinearProgressionParameters true_params(fbp, bgp);
    double log_lik_at_truth = compute_pathway_likelihood(*obs_matrix, row_sum, true_state, true_params);
    cout << "Likelihood at truth: " << log_lik_at_truth << endl;

    // 1. run model selection
    size_t best_l = infer_model_length(random, *obs_matrix, row_sum, max_model_len, num_samples, n_smc_iter, n_mcmc_iter, n_particles, ess, fbp_max, bgp_max, swap_prob, move_type, sim_dat_path, n_threads);

    // 2. run Particle Gibbs: generate parameters and pathways (samples)
    run_particle_gibbs(random, *obs_matrix, row_sum, best_l, n_smc_iter, n_mcmc_iter, n_particles, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_prob, move_type, sim_dat_path, n_threads);
}

void test()
{
    vector<unordered_set<size_t>> a;
    unordered_set<size_t> p1;
    p1.insert(0);
    unordered_set<size_t> p2;
    p2.insert(1);
    a.push_back(p1);
    a.push_back(p2);
    vector<unordered_set<size_t>> b = a;
    a[0].insert(2);
    for (size_t i = 0; i < a.size(); i++)
    {
        cout << i << ": ";
        for (size_t idx : a[i])
        {
            cout << idx << " ";
        }
        cout << endl;
    }
    
    for (size_t i = 0; i < b.size(); i++)
    {
        cout << i << ": ";
        for (size_t idx : a[i])
        {
            cout << idx << " ";
        }
        cout << endl;
    }
}

void test_caching(gsl_rng *random, gsl_matrix &obs, vector<size_t> &row_sum)
{
    size_t n_patients = obs.size1;
    size_t n_genes = obs.size2;
    size_t n_drivers = 3;
    size_t n_pathways = n_drivers + 1;
    LinearProgressionState state(obs, row_sum, n_genes, n_drivers, true);
    
    // randomly assign genes to pathways
    state.sample_from_prior(random);
    
    for (size_t m = 0; m < n_patients; m++) {
        vector<size_t> r(n_pathways);
        fast_matrix_product(obs, m, state.get_pathway_membership(), r);
        
        // compare r vs state's cache version
        const vector<size_t> &cache_r = state.get_cache_at(m);
        
        if (r.size() != cache_r.size()) {
            cerr << "Cache test failed: size mismatch." << endl;
            exit(-1);
        }
        
        for (size_t i = 0; i < r.size(); i++) {
            if (r[i] != cache_r[i]) {
                cerr << "Cache test failed: calculation mismatch." << endl;
                exit(-1);
            }
        }
    }
    cout << "Cache test passed!" << endl;
}

void test_local_caching(gsl_rng *random, gsl_matrix &obs, vector<size_t> row_sum)
{
    size_t n_genes = obs.size2;
    size_t n_drivers = 10;
    LinearProgressionState state(obs, row_sum, n_genes, n_drivers, true);
    
    // randomly assign genes to pathways
    state.sample_from_prior(random);
    
    LinearProgressionParameters params(gsl_ran_flat(random, 0, 0.3), gsl_ran_flat(random, 0, 0.3));

    // compute using cache and compare to without caching
    // test that the same values are attained and less time is required
    
    auto start = std::chrono::high_resolution_clock::now();
    double log_lik_w_cache = compute_pathway_likelihood(obs, row_sum, state, params, true);
    auto end1 = std::chrono::high_resolution_clock::now();
    double log_lik_w_o_cache = compute_pathway_likelihood(obs, row_sum, state, params, false);
    auto end2 = std::chrono::high_resolution_clock::now();
    if (abs(log_lik_w_cache - log_lik_w_o_cache) > 1e-6) {
        cerr << "Error: cache test 2 failed!" << endl;
        exit(-1);
    }
    std::chrono::duration<double> elapsed1 = end1 - start;
    std::chrono::duration<double> elapsed2 = end2 - end1;
    cout << "With cache time elapsed: " << elapsed1.count() << endl;
    cout << "Without cache time elapsed: " << elapsed2.count() << endl;
    cout << "Cache test 2 passed!" << endl;
}

void test_global_caching(gsl_rng *random, gsl_matrix &obs, vector<size_t> row_sum)
{
    size_t n_genes = obs.size2;
    size_t n_drivers = 10;
    size_t n_particles = 100;
    unsigned long seed = gsl_rng_get(random);
    gsl_rng *random1 = generate_random_object(seed);
    gsl_rng *random2 = generate_random_object(seed);
    LinearProgressionModel model(n_genes, n_drivers, 10, 1, obs, row_sum, 0.2, true, MoveType::MH);
    LinearProgressionModelCache model_cache(n_genes, n_drivers, 10, 1, obs, row_sum, 0.2, true, MoveType::MH);
    LinearProgressionParameters params(gsl_ran_flat(random, 0, 0.3), gsl_ran_flat(random, 0, 0.3));
    double logw = 0.0, logw_cache = 0.0;
    std::chrono::duration<double> elapsed_cache1, elapsed_cache2;
    for (size_t i = 0; i < n_particles; i++) {
        auto start1 = std::chrono::high_resolution_clock::now();
        model.propose_initial(random1, logw, params);
        auto end1 = std::chrono::high_resolution_clock::now();
        auto start2 = std::chrono::high_resolution_clock::now();
        model_cache.propose_initial(random2, logw_cache, params);
        auto end2 = std::chrono::high_resolution_clock::now();
        elapsed_cache1 += (end1 - start1);
        elapsed_cache2 += (end2 - start2);
        if (abs(logw - logw_cache) > 1e-6) {
            cerr << "Error: cache version is not correctly implemented." << endl;
            exit(-1);
        }
    }
    cout << "With local caching time elapsed: " << elapsed_cache1.count() << endl;
    cout << "With global caching time elapsed: " << elapsed_cache2.count() << endl;
    cout << "Cache test 3 passed!" << endl;
}

void test_openmp()
{
    omp_set_num_threads(4);
#pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf("%i\n", ID);
    }
}

void test_map()
{
    unordered_map<size_t, unordered_map<size_t, double>> f;
    f[1][2] = 0.1;
    f[2][3] = 0.2;
    f[1][2] = 0.3;
    cout << f[1][2] << ", " << f[2][3] << endl;
}

int main(int argc, char *argv[]) // TODO: take arguments from command line
{
    //test_openmp();
    //test_map();

    unsigned long seed;

    // model length inference parameters
    bool infer_model_len;
    unsigned int max_model_len = 0;
    
    // run pg?
    bool run_pg;
    
    // joint inference of progression pathway and parameters
    unsigned int model_len;

    // smc and pmcmc parameters
    size_t n_mc_samples;
    double ess;
    size_t n_particles;
    size_t n_smc_iter;
    size_t n_mcmc_iter;
    size_t n_pmcmc_iter;
    size_t n_mh_w_gibbs_iter;
    
    // uniform prior on the error probabilities: \delta \sim Uniform(0, fbp_max), \epsilon \sim Uniform(0, bgp_max)
    double fbp_max;
    double bgp_max;
    
    // SMC move performs swap move or
    double swap_move_prob;
    string move_type_str;
    MoveType move_type;
    
    string output_dir;
    
    size_t n_threads;

    // options
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help", "Infer tumor progression pathway.")
    ("data-path,d", po::value<string>(), "Path to observed mutation file (MxK csv file of binary values).")
    ("experiment-path,x", po::value< string >(), "Path to simulated data (for development purposes).")
    ("output_dir,o", po::value<string>( &output_dir )->default_value("."), "Output path. Default will be the current folder.")
    ("seed,s", po::value< unsigned long >( &seed )->default_value(1), "Random seed.")
    ("infer-len,i", po::value< bool >(&infer_model_len)->default_value(false), "Infer best model length.")
    ("run-pg,r", po::value< bool >(&run_pg)->default_value(true), "Run Particle Gibbs to infer pathways. Default is true. Setting -i true and setting -r false if only goal is to infer the model length.")
    ("max-len,L", po::value< unsigned int >( &max_model_len )->default_value(8), "Max model length to use for inferring best model length. If -i option is set to true, this value must be provided.")
    ("model-len,l", po::value< unsigned int >( &model_len )->default_value(5), "Model length to use for inferring Particle Gibbs. If -i option is set to true, this value will inferred.")
    ("n-mc-samples", po::value< size_t >( &n_mc_samples )->default_value(20), "Number of Monte Carlo sampels for model selection.")
    ("ess", po::value< double >( &ess )->default_value(0.5), "A vlue in [0, 1]. ESS for triggering resampling of SMC sampler for model selection.")
    ("n-particles", po::value< size_t >( &n_particles )->default_value(200), "Number of SMC samples.")
    ("n-smc-iter", po::value< size_t >( &n_smc_iter )->default_value(100), "Number of SMC iterations.")
    ("n-mcmc-iter", po::value< size_t >( &n_mcmc_iter )->default_value(5), "Number of MCMC iterations to be used within SMC sampler.")
    ("n-pmcmc-iter", po::value< size_t >( &n_pmcmc_iter )->default_value(1000), "Number of Particle MCMC iterations.")
    ("n-mh-w-gibbs-iter", po::value< size_t >( &n_mh_w_gibbs_iter )->default_value(10), "Number of MH-within-Gibbs iterations for sampling error probabilities.")
    ("bgp-max", po::value< double >( &bgp_max )->default_value(0.3), "A vlue in (0, 1]. Parameterize uniform prior on the parameters: Uniform[0, bgp-max].")
    ("fbp-max", po::value< double >( &fbp_max )->default_value(0.3), "A vlue in (0, 1]. Parameterize uniform prior on the parameters: Uniform[0, fbp-max].")
    ("swap-prob", po::value< double >( &swap_move_prob )->default_value(0.2), "A vlue in [0, 1]. Indiciate probability for pathway swap move probability.")
    ("move-type", po::value< string >( &move_type_str )->default_value("MH"), "One of {MH, GIBBS}. Invalid specification and default is MH.")
    ("n-threads,t", po::value< size_t >( &n_threads )->default_value(4), "Number of threads available.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("seed")) {
        cout << "Random seed: " << vm["seed"].as<size_t>() << endl;
    }
    if (vm.count("infer-model-len")) {
        cout << "Infer model length: " << (infer_model_len ? "true" : "false") << endl;
        if (vm.count("max-len")) {
            cout << "Max model length: " << vm["max-len"].as<unsigned int>() << endl;
        } else {
            cerr << "-i option is set to true. Please provide max model length." << endl;
            exit(-1);
        }
    }
    else {
        cout << "Model length: " << model_len << endl;
    }
    move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;
    
    gsl_rng *random = generate_random_object(seed);

    if (vm.count("experiment")) {
        // development mode -- data file, if provided, will be ignored
        string sim_data_path = vm["experiment"].as<string>();
        run_exp(random, sim_data_path, max_model_len, n_mc_samples, n_smc_iter, n_mcmc_iter, n_particles, ess, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_move_prob, move_type, n_threads);
    } else {
        string dat_file;
        if (vm.count("data-path")) {
            dat_file = vm["data-path"].as<string>();
            cout << "Data file:" << dat_file << ".\n";
        } else {
            cerr << "Data file was not set.\n";
            exit(-1);
        }

        // read the data and compute row sum
        gsl_matrix *obs_matrix = read_data(dat_file, true);
        vector<size_t> row_sum(obs_matrix->size1);
        compute_row_sum(*obs_matrix, row_sum);

        test_caching(random, *obs_matrix, row_sum);
        test_local_caching(random, *obs_matrix, row_sum);
        test_global_caching(random, *obs_matrix, row_sum);
        
        if (infer_model_len) {
            model_len = infer_model_length(random, *obs_matrix, row_sum, max_model_len, n_mc_samples, n_smc_iter, n_mcmc_iter, n_particles, ess, fbp_max, bgp_max, swap_move_prob, move_type, output_dir, n_threads);
            cout << "Inferred model length: " << model_len << endl;
        }

        if (run_pg) {
            cout << "Running Particle Gibbs with model length: " << model_len << endl;
            run_particle_gibbs(random, *obs_matrix, row_sum, model_len, n_smc_iter, n_mcmc_iter, n_particles, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_move_prob, move_type, output_dir, n_threads);
        }
    }

    return 0;
}
