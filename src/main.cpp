//
//  main.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-05-08.
//

#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>
#include "lpm.hpp"

using namespace std;

void draw_stratified_samples(gsl_rng *random, size_t n_mc_samples, double vec[], double max)
{
    double interval = max/n_mc_samples;
    double sum = 0.0;
    for (size_t i = 0; i < n_mc_samples; i++) {
        vec[i] = gsl_ran_flat(random, sum, sum + interval);
        sum += interval;
    }
}

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    /**
     * arguments:
     * mode: model/pg
     * seed
     * model_len: range (i.e., 2-10)
     * n_mc_samples: 30
     * n_pg_iter: 30
     * n_particles: 400
     * n_smc_iter: 100
     * n_kernel_iter: 3
     * n_mh_w_gibbs_iter: 10
     * has_passenger: True
     * swap_prob: 0.1
     * fbp_max: 0
     * bgp_max: 0.2
     * mh_proposal_sd: 0.05
     * n_threads: 8
     * use_lik_tempering: True
     * data_path: data/Experiment2/With_passengers/5/error0.05/rep0/matrix.csv
     **/
    
    string data_path;
    string output_path;
    string mode;
    unsigned long seed;
    unsigned int model_len;
    //unsigned int n_mc_samples;
    unsigned int n_mcmc_iter;
//    unsigned int n_particles;
//    unsigned int n_smc_iter;
//    unsigned int n_kernel_iter;
    unsigned int n_mh_w_gibbs_iter;
    unsigned int thinning_interval;
//    unsigned int n_mc_threads;
//    unsigned int n_smc_threads;
//    bool use_lik_tempering;
    double fbp_max, bgp_max;
    bool has_passenger;
    double swap_prob;
    double mh_proposal_sd;
    double prior_passenger_prob;

    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("data_path,d", po::value<string>(&data_path)->required(), "Specify path to the data.")
    ("output_path,o", po::value<string>(&output_path)->required(), "Specify output path.")
    ("mode,m", po::value<string>(&mode)->required(), "Specify mode: model for model selection and pg for particle Gibbs.")
    ("seed,s", po::value<unsigned long>(&seed)->default_value(1), "Specify random seed.")
    ("model_len,l", po::value<unsigned int>(&model_len)->required(), "Specify model length.")
    //("n_mc_samples,M", po::value<unsigned int>(&n_mc_samples), "Specify number of MC samples to use for model selection.")
    ("n_mcmc_iter,p", po::value<unsigned int>(&n_mcmc_iter), "Specify number of iterations for PG/MCMC.")
    //("n_particles,P", po::value<unsigned int>(&n_particles)->required(), "Specify number of particles to use.")
    //("n_smc_iter,S", po::value<unsigned int>(&n_smc_iter)->required(), "Specify number of SMC iterations.")
    //("n_kernel_iter,k", po::value<unsigned int>(&n_kernel_iter)->default_value(5), "Specify number of iterations of MCMC kernels to apply within a SMC move. Specifying 0 reduces to random walk move.")
    ("n_mh_w_gibbs_iter,G", po::value<unsigned int>(&n_mh_w_gibbs_iter)->default_value(10), "Specify number of iterations of MH within Gibbs for inferring the parameters within PG.")
    //("n_mc_threads,t", po::value<unsigned int>(&n_mc_threads)->default_value(1), "Specify number of threads to use for parallel processing of Monte Carlo samples.")
    //("n_smc_threads", po::value<unsigned int>(&n_smc_threads)->default_value(1), "Specify number of threads to use for SMC algorithm.")
    //("use_lik_tempering,T", po::value<bool>(&use_lik_tempering)->default_value(true), "Specify whether to use likelihood tempering.")
    ("fbp_max,f", po::value<double>(&fbp_max)->default_value(0), "Specify max value for flip-back probability. Specify 0 here to use unified error probability.")
    ("bgp_max,b", po::value<double>(&bgp_max)->default_value(0.2), "Specify max value for background mutation probability. Specifying 0 for fbp_max to interpret this quantity as an error probability.")
    ("has_passenger", po::value<bool>(&has_passenger)->default_value(true), "Specify whether to detect passenger genes. 0: false, 1: true.")
    ("swap_prob", po::value<double>(&swap_prob)->default_value(0.2), "Specify pathway swap probability for SMC move.")
    ("mh_proposal_sd", po::value<double>(&mh_proposal_sd)->default_value(0.05), "Specify standard deviation of Gaussian random walk proposal within MHwGibbs for inferring the parameters (for mode = pg).")
    ("thinning_interval", po::value<unsigned int>(&thinning_interval)->default_value(1), "Specify thinning interval for MCMC.")
    ("prior_passenger_prob", po::value<double>(&prior_passenger_prob)->default_value(0.95), "Specify prior probability of a gene being a passenger.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    
    if (mode == "mcmc") {
        run_mcmc(seed, data_path.c_str(), output_path.c_str(), model_len, n_mcmc_iter, n_mh_w_gibbs_iter, thinning_interval, has_passenger, swap_prob, fbp_max, bgp_max, mh_proposal_sd, prior_passenger_prob);
    } else {
        cerr << "Unknown mode: " << mode << "." << endl;
    }
    
//    if (mode == "pg") {
//        if (!vm.count("n_pg_iter")) {
//            cerr << "Please specify number of PG iterations." << endl;
//        }
//        run_pg(seed, data_path.c_str(), output_path.c_str(), model_len, n_mcmc_iter, n_particles, n_smc_iter, n_kernel_iter, n_mh_w_gibbs_iter, has_passenger, swap_prob, fbp_max, bgp_max, mh_proposal_sd, n_smc_threads);
//    } else if (mode == "model") {
//        if (!vm.count("n_mc_samples")) {
//            cerr << "Please specify number of Monte Carlo samples." << endl;
//        }
//        // sample fbps and bgps via stratified sampling
//        double *fbps = new double[n_mc_samples];
//        double *bgps = new double[n_mc_samples];
//        gsl_rng *random = generate_random_object(seed);
//        draw_stratified_samples(random, n_mc_samples, bgps, bgp_max);
//        if (fbp_max == 0.0) {
//            std::copy(bgps, bgps + n_mc_samples, fbps);
//        }
//        unsigned long new_seed = gsl_rng_get(random);
//        model_selection(new_seed, data_path.c_str(), output_path.c_str(), model_len, n_mc_samples, n_particles, n_smc_iter, n_kernel_iter, has_passenger, swap_prob, fbps, bgps, n_mc_threads, n_smc_threads, use_lik_tempering);
//        gsl_rng_free(random);
//    } else if (mode == "mcmc") {
//        run_mcmc(seed, data_path.c_str(), output_path.c_str(), model_len, n_mcmc_iter, n_mh_w_gibbs_iter, thinning_interval, has_passenger, swap_prob, fbp_max, bgp_max, mh_proposal_sd, prior_passenger_prob);
//    } else {
//        cerr << "Unknown mode: " << mode << "." << endl;
//    }
    return 0;
}
