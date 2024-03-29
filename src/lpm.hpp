//
//  lpm.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-27.
//

#ifndef lpm_h
#define lpm_h

#include <gsl/gsl_matrix.h>

#include "data_util.hpp"
#include "lpm_state.hpp"
#include "lpm_params.hpp"
#include "sampling_utils.hpp"

extern "C" {
    void run_mcmc(long seed,
                  const char *dat_file,
                  const char *output_path,
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
                  double prior_passenger_prob = 0.9);
    
//    void run_pg(long seed,
//                const char *dat_file,
//                const char *output_path,
//                unsigned int model_len,
//                unsigned int n_pg_iter,
//                unsigned int n_particles,
//                unsigned int n_smc_iter,
//                unsigned int n_kernel_iter,
//                unsigned int n_mh_w_gibbs_iter,
//                bool has_passenger,
//                double swap_prob,
//                double fbp_max,
//                double bgp_max,
//                double mh_proposal_sd,
//                bool use_lik_tempering = false,
//                unsigned int n_threads = 4);
//
//    double run_smc(long seed,
//                   const char *dat_file,
//                   unsigned int model_len,
//                   unsigned int n_particles,
//                   unsigned int n_smc_iter,
//                   unsigned int n_kernel_iter,
//                   bool has_passenger,
//                   double swap_prob,
//                   double fbp,
//                   double bgp,
//                   unsigned int *states,
//                   double *log_weights,
//                   bool use_lik_tempering = false,
//                   unsigned int n_threads = 4);

    void process_samples(const string output_path,
                         const gsl_matrix &data,
                         const gsl_matrix *states,
                         vector<double> &posterior_bgps,
                         vector<double> &posterior_fbps,
                         unsigned int n_drivers,
                         bool has_passenger,
                         double prior_passenger_prob);
    void perform_prediction(const char *mcmc_samples_path,
                            const char *data_path,
                            unsigned int model_len,
                            bool has_passenger,
                            unsigned int *true_pathway = 0);
    
    // input fbps and bgps, the parameters to use
    // log_marginal_sum[n]: sum estimation
    // log_marginal_smc[n]: SMC estimation
    // return value: 
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
                           bool use_lik_tempering = false);

    double compute_likelihood(const char *dat_file,
                              unsigned int *pathway,
                              unsigned int model_len,
                              unsigned int n_genes,
                              bool has_passenger,
                              double fbp,
                              double bgp);

    void generate_data(long seed,
                       const char *output_path,
                       unsigned int n_patients,
                       unsigned int model_len,
                       unsigned int n_genes,
                       double fbp,
                       double bgp);
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
                          double prior_passenger_prob = 0.9);

//ParticleGibbs<LinearProgressionState, LinearProgressionParameters> run_pg_from_matrix(
//                        long seed,
//                        gsl_matrix *obs_matrix,
//                        unsigned int model_len,
//                        unsigned int n_pg_iter,
//                        unsigned int n_particles,
//                        unsigned int n_smc_iter,
//                        unsigned int n_kernel_iter,
//                        unsigned int n_mh_w_gibbs_iter,
//                        bool has_passenger,
//                        double swap_prob,
//                        double fbp_max,
//                        double bgp_max,
//                        vector<double> *fbps = 0,
//                        vector<double> *bgps = 0,
//                        double mh_proposal_sd = 0.05,
//                        bool use_lik_tempering = false,
//                        const char *output_path = nullptr,
//                        unsigned int n_threads = 4);
//
//double run_smc_from_matrix(long seed,
//                           gsl_matrix *data_matrix,
//                           unsigned int model_len,
//                           unsigned int n_particles,
//                           unsigned int n_smc_iter,
//                           unsigned int n_kernel_iter,
//                           bool has_passenger,
//                           double swap_prob,
//                           double fbp,
//                           double bgp,
//                           unsigned int *states,
//                           double *log_weights,
//                           bool use_lik_tempering = false,
//                           unsigned int n_threads = 4);

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
                       bool use_lik_tempering = false);

double compute_likelihood_from_matrix(const gsl_matrix *obs_matrix,
                                      unsigned int *pathway,
                                      unsigned int model_len,
                                      unsigned int n_genes,
                                      bool has_passenger,
                                      double fbp,
                                      double bgp);

vector<double> compute_likelihood_per_row(const gsl_matrix *obs_matrix,
                                      unsigned int *pathway,
                                      unsigned int model_len,
                                      unsigned int n_genes,
                                      bool has_passenger,
                                      double fbp,
                                      double bgp,
                                      unsigned int stage);

#endif /* lpm_h */
