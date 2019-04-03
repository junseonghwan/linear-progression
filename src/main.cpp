//#include <chrono>
//#include <iostream>
//#include <math.h>
//#include <stdio.h>
//#include <thread>
//#include <unordered_map>
//#include <vector>
//#include <string>
//
//#include <omp.h>
//
//#include <boost/filesystem.hpp>
//#include <boost/program_options.hpp>
//
//#include <gsl/gsl_statistics.h>
//
//#include <spf/particle_population.hpp>
//#include <spf/pg.hpp>
//#include <spf/sampling_utils.hpp>
//#include <spf/csmc.hpp>
//#include <spf/smc.hpp>
//#include <spf/spf.hpp>
//
//#include "data_util.hpp"
//#include "lpm_likelihood.hpp"
//#include "lpm_model.hpp"
//#include "lpm_params.hpp"
//#include "lpm_pg_proposal.hpp"
//#include "lpm_state.hpp"
//
//using namespace std;

//
////void run_particle_gibbs(gsl_rng *random, gsl_matrix &obs_matrix, vector<size_t> &row_sum,
////                        size_t target_model_length,
////                        size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles,
////                        size_t n_pmcmc_iter, size_t n_mh_gibbs_iter,
////                        double fbp_max, double bgp_max,
////                        double swap_prob, MoveType move_type,
////                        string output_dir, size_t n_threads,
////                        bool has_passenger)
////{
////    size_t n_patients = obs_matrix.size1;
////    size_t n_genes = obs_matrix.size2;
////
////    if (target_model_length >= n_genes)
////    {
////        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
////        exit(-1);
////    }
////
////    // smc settings
////    SMCOptions smc_options;
////    smc_options.num_particles = n_particles;
////    smc_options.ess_threshold = 1.0;
////    smc_options.resample_last_round = false;
////
////    // pmcmc options
////    size_t random_seed = gsl_rng_get(random);
////
////    shared_ptr<LinearProgressionState> prev_state_ptr;
////    for (size_t l = 1; l <= target_model_length; l++) {
////        cout << "Running with model length " << l << endl;
////        LinearProgressionModel smc_model(n_genes, l, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum, swap_prob, has_passenger, move_type);
////        if (l > 1) {
////            smc_model.set_initial_state(random, prev_state_ptr.get());
////        }
////
////        smc_options.main_seed = gsl_rng_get(random);
////        smc_options.resampling_seed = gsl_rng_get(random);
////
////        ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
////        // for l < target, we don't need to generate too many samples as convergence is not the concern
////        PMCMCOptions pmcmc_options(random_seed, round(n_pmcmc_iter));
////        LPMParamProposal pg_proposal(n_patients, n_mh_gibbs_iter, fbp_max, bgp_max);
////        ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(pmcmc_options, csmc, pg_proposal);
////
////        auto start = std::chrono::high_resolution_clock::now();
////        pg.run();
////        auto end = std::chrono::high_resolution_clock::now();
////        std::chrono::duration<double> elapsed = end - start;
////        cout << elapsed.count() << " seconds elapsed." <<  endl;
////
////        // get the last genealogy, last state to initialize the next step
////        vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states = pg.get_states();
////        ParticleGenealogy<LinearProgressionState> *genealogy = states.at(states.size() - 1).get();
////        const LinearProgressionState &state = genealogy->get_state_at(genealogy->size() - 1);
////        prev_state_ptr.reset(new LinearProgressionState(state)); // make a copy
////        cout << "Current solution:\n" << state.to_string() << endl;
////
////        if (l == target_model_length) {
////            vector<shared_ptr<LinearProgressionParameters> > &params = pg.get_parameters();
////            write_pg_output(output_dir, states, params);
////        }
////    }
////    // output parameters and states for the targe model length
////
////}
//
//void run_particle_gibbs(gsl_rng *random, gsl_matrix &obs_matrix, vector<size_t> &row_sum,
//                        size_t target_model_length,
//                        size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles,
//                        size_t n_pmcmc_iter, size_t n_mh_gibbs_iter,
//                        double fbp_max, double bgp_max,
//                        double swap_prob, MoveType move_type,
//                        string output_dir, size_t n_threads,
//                        bool has_passenger)
//{
//    size_t n_patients = obs_matrix.size1;
//    size_t n_genes = obs_matrix.size2;
//
//    if (target_model_length >= n_genes)
//    {
//        cerr << "Maximum model length cannot be larger than the number of genes." << endl;
//        exit(-1);
//    }
//
//    cout << "Running with model length " << target_model_length << endl;
//
//    // smc settings
//    SMCOptions smc_options;
//    smc_options.num_particles = n_particles;
//    smc_options.ess_threshold = 1.0;
//    smc_options.resample_last_round = false;
//    smc_options.main_seed = gsl_rng_get(random);
//    smc_options.resampling_seed = gsl_rng_get(random);
//
//    // pmcmc options
//    size_t random_seed = gsl_rng_get(random);
//
//    LinearProgressionModel smc_model(n_genes, target_model_length, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum, swap_prob, has_passenger, move_type);
//
//    
//    ConditionalSMC<LinearProgressionState, LinearProgressionParameters> csmc(smc_model, smc_options);
//    PMCMCOptions pmcmc_options(random_seed, round(n_pmcmc_iter));
//    LPMParamProposal pg_proposal(n_patients, n_mh_gibbs_iter, fbp_max, bgp_max);
//    ParticleGibbs<LinearProgressionState, LinearProgressionParameters> pg(pmcmc_options, csmc, pg_proposal);
//    
//    auto start = std::chrono::high_resolution_clock::now();
//    pg.run();
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    cout << elapsed.count() << " seconds elapsed." <<  endl;
//
//    // get the last genealogy, last state to initialize the next step
//    vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states = pg.get_states();
//    ParticleGenealogy<LinearProgressionState> *genealogy = states.at(states.size() - 1).get();
//    const LinearProgressionState &state = genealogy->get_state_at(genealogy->size() - 1);
//    cout << "Current solution:\n" << state.to_string() << endl;
//
//    //vector<shared_ptr<LinearProgressionParameters> > &params = pg.get_parameters();
//    write_pg_output(output_dir, states, pg_proposal);
//}
//
//
//size_t infer_model_length(gsl_rng *random, gsl_matrix &obs_matrix, vector<size_t> &row_sum,
//                            size_t max_model_len, size_t num_samples,
//                            size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles, size_t ess,
//                            double fbp_max, double bgp_max,
//                            double swap_prob, MoveType move_type,
//                            string output_dir, size_t num_threads,
//                            bool has_passenger)
//{
//    size_t n_patients = obs_matrix.size1;
//    size_t n_genes = obs_matrix.size2;
//    
//    cout << "Inferring model length: " << endl;
//    cout << "# Patients: " << n_patients << endl;
//    cout << "# Genes: " << n_genes << endl;
//    cout << "# Particles: " << n_particles << endl;
//    cout << "# MC Samples: " << num_samples << endl;
//    cout << "# SMC iters: " << n_smc_iter << endl;
//    cout << "# MCMC iters: " << n_mcmc_iter << endl;
//    cout << "Allocate passenger pathway: " << has_passenger << endl;
//
//    if (max_model_len >= n_genes)
//    {
//        cerr << "Maximum model length must be smaller than the number of genes." << endl;
//        exit(-1);
//    }
//
//    //double log_marginals[max_model_len][num_samples];
//    vector<vector<double>> log_marginals;
//    vector<double> model_evidences;
//    double best_marginal_log_lik = DOUBLE_NEG_INF;
//    size_t best_l = 0;
//
//    // sample the parameters in advance and use the same values throughout
//    vector<double> bgps;
//    vector<double> fbps;
//    for (size_t i = 0; i < num_samples; i++) {
//        double bgp = gsl_ran_flat(random, 0.0, bgp_max);
//        double fbp;
//        if (fbp_max > 0.0) {
//            fbp = gsl_ran_flat(random, 0.0, fbp_max);
//        } else {
//            fbp = bgp; // unify the error probs
//        }
//        bgps.push_back(bgp);
//        fbps.push_back(fbp);
//    }
//
//    for (size_t l = 2; l <= max_model_len; l++) {
//        cout << "Model length: " << l << endl;
//        vector<double> log_marginal_vector(num_samples);
//        
//        auto start = std::chrono::high_resolution_clock::now();
//
//#pragma omp parallel for
//    {
//        for (size_t i = 0; i < num_samples; i++) {
//            LinearProgressionModel model(n_genes, l, n_smc_iter, n_mcmc_iter, obs_matrix, row_sum, swap_prob, has_passenger, move_type);
//            LinearProgressionParameters params(fbps[i], bgps[i]);
//            SMCOptions smc_options;
//            smc_options.ess_threshold = ess;
//            smc_options.num_particles = n_particles;
//            smc_options.resample_last_round = false;
//            smc_options.debug = false;
//            smc_options.num_threads = 1;
//
//            ConditionalSMC<LinearProgressionState, LinearProgressionParameters> smc(model, smc_options);
//            smc.initialize(params);
//            log_marginal_vector[i] = smc.get_log_marginal_likelihood();
//        }
//    }
//        
//        auto end = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double> elapsed = end - start;
//        cout << "Total running time: "  << elapsed.count() << endl;
//        
//        double model_evidence = DOUBLE_NEG_INF;
//        for (size_t i = 0; i < num_samples; i++) {
//            model_evidence = log_add(model_evidence, log_marginal_vector[i]);
//        }
//        model_evidence -= log(num_samples);
//        cout << "Model evidence: " << model_evidence << endl;
//        if (model_evidence > best_marginal_log_lik) {
//            best_marginal_log_lik = model_evidence;
//            best_l = l;
//        }
//        
//        log_marginals.push_back(log_marginal_vector);
//        model_evidences.push_back(model_evidence);
//
//    }
//    
//    write_model_selection_output_to_file(output_dir, n_patients, model_evidences, fbps, bgps, log_marginals);
//
//    return best_l;
//}
//
//void run_exp(gsl_rng *random, string exp_path,
//             size_t max_model_len, size_t num_samples,
//             size_t n_smc_iter, size_t n_mcmc_iter, size_t n_particles, size_t ess,
//             size_t n_pmcmc_iter, size_t n_mh_w_gibbs_iter,
//             double fbp_max, double bgp_max,
//             double swap_prob, MoveType move_type,
//             size_t n_threads, bool has_passenger)
//{
//    // for each folder rep[0-9]+:
//    // 1. infer the model length for patients {20, 50, 100, 200, 500, 1000}.
//    // 2. run PG with inferred model length.
//    // 3. run PG with true model length: can find it from generative_mem_mat.csv.
//    
//    vector<size_t> num_patients_vector;
////    num_patients_vector.push_back(20);
////    num_patients_vector.push_back(50);
////    num_patients_vector.push_back(100);
////    num_patients_vector.push_back(200);
////    num_patients_vector.push_back(500);
//    num_patients_vector.push_back(1000);
//    
//    size_t rep_idx = 20;
//    while (true) {
//        string sim_dat_path = exp_path + "/rep" + to_string(rep_idx++);
//
//        // check if the directory exists
//        if (!path_exists(sim_dat_path))
//            break;
//
//        // path to files: data
//        // ground truth: pathway, patient stages, params, model length used for generating the data
//        string obs_file = sim_dat_path + "/matrix.csv";
//        string ground_truth_file = sim_dat_path + "/generative_mem_mat.csv";
//        string stages_file = sim_dat_path + "/patients_stages.csv";
//        string params_file = sim_dat_path + "/parameters.csv";
//
//        // read the data
//        gsl_matrix *obs_matrix = read_data(obs_file);
//        size_t n_patients = obs_matrix->size1;
//        size_t n_genes = obs_matrix->size2;
//
//        // compute the row sum
//        vector<size_t> row_sum(n_patients);
//        compute_row_sum(*obs_matrix, row_sum);
//
//        // read ground truth model length
//        gsl_matrix *pathway_membership = read_data(ground_truth_file);
//        size_t true_model_length = pathway_membership->size2;
//
//        // read the ground truth pathway
//        LinearProgressionState true_state(*obs_matrix, row_sum, n_genes, true_model_length, has_passenger);
//        for (size_t g = 0; g < pathway_membership->size1; g++) {
//            size_t pathway = 0;
//            for (size_t k = 0; k < pathway_membership->size2; k++) {
//                if (gsl_matrix_get(pathway_membership, g, k) == 1) {
//                    pathway = k;
//                    break;
//                }
//            }
//            true_state.update_pathway_membership(g, pathway);
//        }
//
//        // read the parameters used for generating the data
//        double fbp = 0.0;
//        double bgp = 0.0;
//        read_error_params(params_file, fbp, bgp);
//
//        // extract the first N lines of obs_matrix and row_sum
//        
//        // compute likelihood of the truth
//        LinearProgressionParameters true_params(fbp, bgp);
//        double log_lik_at_truth = compute_pathway_likelihood(*obs_matrix, row_sum, true_state, true_params);
//        cout << "Likelihood at truth: " << log_lik_at_truth << endl;
//        
//        // testing: read in the gold standard likelihood and compare
//        string gold_standard_lik_file = sim_dat_path + "/gold_standard.csv";
//        gsl_matrix *lik_gold_standard = read_csv(gold_standard_lik_file, false);
//        if (obs_matrix->size1 == gsl_matrix_get(lik_gold_standard, lik_gold_standard->size1 - 1, 0)) {
//            double gold_lik = gsl_matrix_get(lik_gold_standard, lik_gold_standard->size1 - 1, 1);
//            if (abs(gold_lik - log_lik_at_truth) > 1e-3) {
//                cerr << "Error: likelihood computation differs from the gold standard." << endl;
//                exit(-1);
//            }
//        }
//
//        // create directory
//        boost::filesystem::path output_path = sim_dat_path + "/inference";
//        boost::filesystem::create_directory(output_path);
//        for (size_t i = 0; i < num_patients_vector.size(); i++)
//        {
//            gsl_matrix *data_matrix = get_first_n_lines(obs_matrix, num_patients_vector[i]);
//            // 1. run model selection
//            size_t best_l = infer_model_length(random, *data_matrix, row_sum, max_model_len, num_samples, n_smc_iter, n_mcmc_iter, n_particles, ess, fbp_max, bgp_max, swap_prob, move_type, output_path.string(), n_threads, has_passenger);
//
//            // 2. run PG using true model length
////            boost::filesystem::path pg_output_path = sim_dat_path + "/inference/npatients" + to_string(num_patients_vector[i]);
////            boost::filesystem::create_directory(pg_output_path);
////            run_particle_gibbs(random, *data_matrix, row_sum, true_model_length, n_smc_iter, n_mcmc_iter, n_particles, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_prob, move_type, pg_output_path.string(), n_threads, has_passenger);
////
////            // 3. run PG using inferred model length if different from the truth
////            if (best_l != true_model_length) {
////                pg_output_path = sim_dat_path + "/inference/pg" + to_string(best_l);
////                boost::filesystem::create_directory(pg_output_path);
////                run_particle_gibbs(random, *data_matrix    , row_sum, best_l, n_smc_iter, n_mcmc_iter, n_particles, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_prob, move_type, pg_output_path.string(), n_threads, has_passenger);
////            }
//
//            delete data_matrix;
//        }
//    }
//}
//
//void test()
//{
//    vector<unordered_set<size_t>> a;
//    unordered_set<size_t> p1;
//    p1.insert(0);
//    unordered_set<size_t> p2;
//    p2.insert(1);
//    a.push_back(p1);
//    a.push_back(p2);
//    vector<unordered_set<size_t>> b = a;
//    a[0].insert(2);
//    for (size_t i = 0; i < a.size(); i++)
//    {
//        cout << i << ": ";
//        for (size_t idx : a[i])
//        {
//            cout << idx << " ";
//        }
//        cout << endl;
//    }
//
//    for (size_t i = 0; i < b.size(); i++)
//    {
//        cout << i << ": ";
//        for (size_t idx : a[i])
//        {
//            cout << idx << " ";
//        }
//        cout << endl;
//    }
//}
//
//void test_caching(gsl_rng *random, gsl_matrix &obs, vector<size_t> &row_sum)
//{
//    size_t n_patients = obs.size1;
//    size_t n_genes = obs.size2;
//    size_t n_drivers = 3;
//    size_t n_pathways = n_drivers + 1;
//    LinearProgressionState state(obs, row_sum, n_genes, n_drivers, true);
//    
//    // randomly assign genes to pathways
//    state.sample_from_prior(random);
//    
//    for (size_t m = 0; m < n_patients; m++) {
//        vector<size_t> r(n_pathways);
//        fast_matrix_product(obs, m, state.get_pathway_membership(), r);
//        
//        // compare r vs state's cache version
//        const vector<size_t> &cache_r = state.get_cache_at(m);
//        
//        if (r.size() != cache_r.size()) {
//            cerr << "Cache test failed: size mismatch." << endl;
//            exit(-1);
//        }
//        
//        for (size_t i = 0; i < r.size(); i++) {
//            if (r[i] != cache_r[i]) {
//                cerr << "Cache test failed: calculation mismatch." << endl;
//                exit(-1);
//            }
//        }
//    }
//    cout << "Cache test passed!" << endl;
//}
//
////void test_local_caching(gsl_rng *random, gsl_matrix &obs, vector<size_t> row_sum)
////{
////    size_t n_genes = obs.size2;
////    size_t n_drivers = 10;
////    LinearProgressionState state(obs, row_sum, n_genes, n_drivers, true);
////
////    // randomly assign genes to pathways
////    state.sample_from_prior(random);
////
////    LinearProgressionParameters params(gsl_ran_flat(random, 0, 0.3), gsl_ran_flat(random, 0, 0.3));
////
////    // compute using cache and compare to without caching
////    // test that the same values are attained and less time is required
////
////    auto start = std::chrono::high_resolution_clock::now();
////    double log_lik_w_cache = compute_pathway_likelihood(obs, row_sum, state, params, true);
////    auto end1 = std::chrono::high_resolution_clock::now();
////    double log_lik_w_o_cache = compute_pathway_likelihood(obs, row_sum, state, params, false);
////    auto end2 = std::chrono::high_resolution_clock::now();
////    if (abs(log_lik_w_cache - log_lik_w_o_cache) > 1e-6) {
////        cerr << "Error: cache test 2 failed!" << endl;
////        exit(-1);
////    }
////    std::chrono::duration<double> elapsed1 = end1 - start;
////    std::chrono::duration<double> elapsed2 = end2 - end1;
////    cout << "With cache time elapsed: " << elapsed1.count() << endl;
////    cout << "Without cache time elapsed: " << elapsed2.count() << endl;
////    cout << "Cache test 2 passed!" << endl;
////}
////
////void test_global_caching(gsl_rng *random, gsl_matrix &obs, vector<size_t> row_sum)
////{
////    size_t n_genes = obs.size2;
////    size_t n_drivers = 10;
////    size_t n_particles = 100;
////    unsigned long seed = gsl_rng_get(random);
////    gsl_rng *random1 = generate_random_object(seed);
////    gsl_rng *random2 = generate_random_object(seed);
////    LinearProgressionModel model(n_genes, n_drivers, 10, 1, obs, row_sum, 0.2, true, MoveType::MH);
////    LinearProgressionModelCache model_cache(n_genes, n_drivers, 10, 1, obs, row_sum, 0.2, true, MoveType::MH);
////    LinearProgressionParameters params(gsl_ran_flat(random, 0, 0.3), gsl_ran_flat(random, 0, 0.3));
////    double logw = 0.0, logw_cache = 0.0;
////    std::chrono::duration<double> elapsed_cache1, elapsed_cache2;
////    for (size_t i = 0; i < n_particles; i++) {
////        auto start1 = std::chrono::high_resolution_clock::now();
////        model.propose_initial(random1, logw, params);
////        auto end1 = std::chrono::high_resolution_clock::now();
////        auto start2 = std::chrono::high_resolution_clock::now();
////        model_cache.propose_initial(random2, logw_cache, params);
////        auto end2 = std::chrono::high_resolution_clock::now();
////        elapsed_cache1 += (end1 - start1);
////        elapsed_cache2 += (end2 - start2);
////        if (abs(logw - logw_cache) > 1e-6) {
////            cerr << "Error: cache version is not correctly implemented." << endl;
////            exit(-1);
////        }
////    }
////    cout << "With local caching time elapsed: " << elapsed_cache1.count() << endl;
////    cout << "With global caching time elapsed: " << elapsed_cache2.count() << endl;
////    cout << "Cache test 3 passed!" << endl;
////}
//
//void test_openmp()
//{
//    omp_set_num_threads(4);
//#pragma omp parallel
//    {
//        int ID = omp_get_thread_num();
//        printf("%i\n", ID);
//    }
//}
//
//void test_map()
//{
//    unordered_map<size_t, unordered_map<size_t, double>> f;
//    f[1][2] = 0.1;
//    f[2][3] = 0.2;
//    f[1][2] = 0.3;
//    cout << f[1][2] << ", " << f[2][3] << endl;
//}
//
//int main(int argc, char *argv[]) // TODO: take arguments from command line
//{
//    //test_openmp();
//    //test_map();
//
//    unsigned long seed;
//
//    // model length inference parameters
//    bool infer_model_len;
//    unsigned int max_model_len = 0;
//    
//    // run pg?
//    bool run_pg;
//    
//    // joint inference of progression pathway and parameters
//    unsigned int model_len;
//
//    // smc and pmcmc parameters
//    size_t n_mc_samples;
//    double ess;
//    size_t n_particles;
//    size_t n_smc_iter;
//    size_t n_mcmc_iter;
//    size_t n_pmcmc_iter;
//    size_t n_mh_w_gibbs_iter;
//    
//    // uniform prior on the error probabilities: \delta \sim Uniform(0, fbp_max), \epsilon \sim Uniform(0, bgp_max)
//    double fbp_max;
//    double bgp_max;
//    
//    // SMC move performs swap move or
//    double swap_move_prob;
//    bool has_passenger;
//    string move_type_str;
//    MoveType move_type;
//
//    string output_dir;
//
//    size_t n_threads;
//
//    // options
//    namespace po = boost::program_options;
//    po::options_description desc("Allowed options");
//    desc.add_options()
//    ("help", "Infer tumor progression pathway.")
//    ("bgp-max,b", po::value< double >( &bgp_max )->default_value(0.1), "A value in (0, 1]. Parameterize uniform prior on the parameters: Uniform[0, bgp-max]. Default 0.3.")
//    ("data-path,d", po::value<string>(), "Path to observed mutation file. MxK matrix with 0's and 1's delineated by comma.")
//    ("ess,e", po::value< double >( &ess )->default_value(1.0), "A vlue in [0, 1]. ESS for triggering resampling of SMC sampler for model selection.")
//    ("fbp-max,f", po::value< double >( &fbp_max )->default_value(0), "A value in (0, 1]. Parameterize uniform prior on the parameters: Uniform[0, fbp-max]. Default 0.")
//    ("n-pmcmc-iter,g", po::value< size_t >( &n_pmcmc_iter )->default_value(10), "Number of Particle Gibbs iterations. Default 10.")
//    ("infer-len,i", po::value< bool >(&infer_model_len)->default_value(false), "Infer best model length. Default false.")
//    ("n-mcmc-iter,k", po::value< size_t >( &n_mcmc_iter )->default_value(1), "Number of MCMC kernels to apply within SMC sampler.")
//    ("model-len,l", po::value< unsigned int >( &model_len )->default_value(5), "Model length to use for inferring Particle Gibbs. If -i true, then this value will be ignored. Default 5.")
//    ("max-len,L", po::value< unsigned int >( &max_model_len )->default_value(8), "Max model length to use for inferring best model length. Default value is 8.")
//    ("n-mc-samples,m", po::value< size_t >( &n_mc_samples )->default_value(20), "Number of Monte Carlo sampels for model selection. Default 20.")
//    ("move-type,M", po::value< string >( &move_type_str )->default_value("MH"), "One of {MH, GIBBS}. Invalid specification and default is MH. Default MH.")
//    ("n-particles,n", po::value< size_t >( &n_particles )->default_value(100), "Number of SMC samples. Default 100.")
//    ("output_dir,o", po::value<string>( &output_dir )->default_value("."), "Output path. Default will be the current folder.")
//    ("n-mh-w-gibbs-iter,p", po::value< size_t >( &n_mh_w_gibbs_iter )->default_value(10), "Number of MH-within-Gibbs iterations for sampling error probabilities. Default 10.")
//    ("run-pg,r", po::value< bool >(&run_pg)->default_value(true), "Run Particle Gibbs to infer pathways. Ex: -i true -r false will infer the model length but not run PG. If -l is provided, it will use that value. Default true.")
//    ("seed,s", po::value< unsigned long >( &seed )->default_value(1), "Random seed.")
//    ("n-smc-iter,S", po::value< size_t >( &n_smc_iter ), "Number of SMC iterations. Default will be number of genes x 2.")
//    ("n-threads,t", po::value< size_t >( &n_threads )->default_value(4), "Number of threads available. Default 4.")
//    ("swap-prob,w", po::value< double >( &swap_move_prob )->default_value(0.2), "A value in [0, 1]. Indiciate probability for pathway swap move probability. Default 0.2.")
//    ("has_passenger,H", po::value< bool >( &has_passenger )->default_value(true), "A boolean. Indiciates whether to allocate passenger pathway. Default: true")
//    ("experiment-path,x", po::value< string >(), "Path to simulated data (for development purposes). Should contain directories labelled rep[0-9]+.")
//    ;
//
//    po::variables_map vm;
//    po::store(po::parse_command_line(argc, argv, desc), vm);
//    po::notify(vm);
//
//    if (vm.count("help")) {
//        cout << desc << "\n";
//        return 1;
//    }
//
//    if (vm.count("seed")) {
//        cout << "Random seed: " << vm["seed"].as<size_t>() << endl;
//    }
//    move_type = move_type_str == "GIBBS" ? MoveType::GIBBS : MoveType::MH;
//
//    gsl_rng *random = generate_random_object(seed);
//    
//    if (vm.count("experiment-path")) {
//        if (!vm.count("n-smc-iter")) {
//            cerr << "# of SMC iterations is not specified." << endl;
//            exit(-1);
//        }
//
//        // development mode -- data file, if provided, will be ignored
//        string sim_data_path = vm["experiment-path"].as<string>();
//        run_exp(random, sim_data_path, max_model_len, n_mc_samples, n_smc_iter, n_mcmc_iter, n_particles, ess, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_move_prob, move_type, n_threads, has_passenger);
//    } else {
//        string dat_file;
//        if (vm.count("data-path")) {
//            dat_file = vm["data-path"].as<string>();
//            cout << "Data file:" << dat_file << ".\n";
//        } else {
//            cerr << "Data file was not set.\n";
//            exit(-1);
//        }
//        
//        if (vm.count("infer-model-len")) {
//            cout << "Infer model length: " << (infer_model_len ? "true" : "false") << endl;
//            if (vm.count("max-len")) {
//                cout << "Max model length: " << vm["max-len"].as<unsigned int>() << endl;
//            } else {
//                cerr << "-i option is set to true. Please provide max model length." << endl;
//                exit(-1);
//            }
//        }
//        else {
//            cout << "Model length: " << model_len << endl;
//        }
//
//        // read the data and compute row sum
//        gsl_matrix *obs_matrix = read_data(dat_file, true);
//        vector<size_t> row_sum(obs_matrix->size1);
//        compute_row_sum(*obs_matrix, row_sum);
//
//        //test_caching(random, *obs_matrix, row_sum);
//        //test_local_caching(random, *obs_matrix, row_sum);
//        //test_global_caching(random, *obs_matrix, row_sum);
//
//        size_t n_genes = obs_matrix->size2;
//        if (!vm.count("n-smc-iter")) {
//            n_smc_iter = n_genes * 2;
//        }
//        
//        if (infer_model_len) {
//            model_len = infer_model_length(random, *obs_matrix, row_sum, max_model_len, n_mc_samples, n_smc_iter, n_mcmc_iter, n_particles, ess, fbp_max, bgp_max, swap_move_prob, move_type, output_dir, n_threads, has_passenger);
//            cout << "Inferred model length: " << model_len << endl;
//        }
//
//        if (run_pg) {
//            cout << "Running Particle Gibbs with model length: " << model_len << endl;
//            run_particle_gibbs(random, *obs_matrix, row_sum, model_len, n_smc_iter, n_mcmc_iter, n_particles, n_pmcmc_iter, n_mh_w_gibbs_iter, fbp_max, bgp_max, swap_move_prob, move_type, output_dir, n_threads, has_passenger);
//        }
//    }
//
//    return 0;
//}

#include "lpm.hpp"

#include <iostream>
#include <stdio.h>
#include <string>

#include <boost/filesystem.hpp>

#include <gsl/gsl_randist.h>

#include <spf/sampling_utils.hpp>

#include "data_util.hpp"

using namespace std;

//int main()
//{
//    const char* output_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/raphael/25genes/error0.05/rep1/pg/";
//    const char* input_path = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/raphael/25genes/error0.05/rep1/matrix.csv";
//    
//    boost::filesystem::path dir(output_path);
//    if (!boost::filesystem::exists(dir)) {
//        if (boost::filesystem::create_directory(dir))
//            std::cout << dir << "....Successfully Created !" << std::endl;
//    }
//    
////    size_t n_particles = 100;
////    size_t n_genes = 25;
////    size_t model_len = 5;
//    
////    run_pg(1, input_path, output_path, 5, 3, n_particles, 25, 1, 20, false, 0.2, 0, 0.2);
////    unsigned int *states = new unsigned int[n_particles*n_genes];
////    double *log_weights = new double[n_particles];
////    run_smc(1, input_path, 5, n_particles, 25, 1, false, 0.2, 0.01, 0.01, states, log_weights);
////    for (size_t i = 0; i < n_particles; i++) {
////        cout << "(";
////        for (size_t j = 0; j < n_genes; j++) {
////            cout << states[i*n_genes + j] << " ";
////        }
////        cout << "), " << log_weights[i] <<  endl;
////    }
//
////    size_t n_mc_samples = 50;
////    size_t n_threads = 8;
////    double *fbps = new double[n_mc_samples];
////    double *bgps = new double[n_mc_samples];
////    gsl_rng *random = generate_random_object(132);
////    for (size_t i = 0; i < n_mc_samples; i++) {
////        bgps[i] = gsl_ran_flat(random, 0.0, 0.1);
////        fbps[i] = bgps[i];
////    }
////    //double ret = model_selection(1, input_path, 0, model_len, n_mc_samples, n_particles, n_genes, 1, false, 0.2, fbps, bgps, n_threads);
////    double *log_marginals = new double[n_mc_samples];
////    double ret = model_selection(1, input_path, model_len, n_mc_samples, n_particles, n_genes, 1, false, 0.2, fbps, bgps, n_threads, log_marginals);
////    cout << ret << endl;
////    for (size_t i = 0; i < n_mc_samples; i++) {
////        cout << bgps[i] << ", " << log_marginals[i] << endl;
////    }
//
//    // path to files: data
//    // ground truth: pathway, patient stages, params, model length used for generating the data
//    string ground_truth_file = "/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/linear-progression/data/raphael/25genes/error0.05/rep1/generative_mem_mat.csv";
//    
//    // read the data
//    gsl_matrix *obs_matrix = read_data(input_path, false);
//    size_t n_patients = obs_matrix->size1;
//    size_t n_genes = obs_matrix->size2;
//    
//    // compute the row sum
//    vector<size_t> row_sum(n_patients);
//    compute_row_sum(*obs_matrix, row_sum);
//    
//    // read ground truth model length
//    gsl_matrix *pathway_membership_matrix = read_data(ground_truth_file);
//    size_t true_model_length = pathway_membership_matrix->size2;
//    
//    // read the ground truth pathway
//    unsigned int *pathway_membership = new unsigned int[n_genes];
//    for (size_t g = 0; g < pathway_membership_matrix->size1; g++) {
//        for (size_t k = 0; k < pathway_membership_matrix->size2; k++) {
//            if (gsl_matrix_get(pathway_membership_matrix, g, k) == 1) {
//                pathway_membership[g] = k;
//                break;
//            }
//        }
//    }
//    
//    double log_lik = compute_likelihood(input_path, pathway_membership, true_model_length, n_genes, false, 0.05, 0.05);
//    cout << "log lik: " << log_lik << endl;
//
//    return 0;
//}

