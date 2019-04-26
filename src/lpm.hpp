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

extern "C" {
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
                double mh_proposal_sd);

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
                   bool is_lik_tempered = false);

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
                           double *ret_smc);
    
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
                        vector<double> *fbps = 0,
                        vector<double> *bgps = 0,
                        double mh_proposal_sd = 0.05,
                        const char *output_path = nullptr);

double run_smc_from_matrix(long seed,
                           gsl_matrix *data_matrix,
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
                           bool is_lik_tempered = false);

double compute_likelihood_from_matrix(gsl_matrix *obs_matrix,
                                      unsigned int *pathway,
                                      unsigned int model_len,
                                      unsigned int n_genes,
                                      bool has_passenger,
                                      double fbp,
                                      double bgp);

#endif /* lpm_h */
