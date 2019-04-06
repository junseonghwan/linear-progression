//
//  lpm.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-27.
//

#ifndef lpm_h
#define lpm_h

#include <gsl/gsl_matrix.h>

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
                double bgp_max);

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
                    double *log_weights);

    
    // let k = model_len
    // let \theta_i = (fbp_i, bgp_i) be a sample from prior p(\theta)
    // run SMC to generate pathway samples {X_j} given \theta_i and k
    // compute: \sum_i (\sum_j p(y | X_j, \theta_i, k) p(X_j | k))
    // this is to approximate f(k) = \int_{\theta} \sum_x  p(y | x, \theta, k) p(x | k) p(\theta) d\theta
    // returns log(\hat{f}(k))
    double model_selection(long seed,
                           const char *dat_file,
                           //const char *output_file,
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
    
    double compute_likelihood_from_matrix(gsl_matrix *obs_matrix,
                                          unsigned int *pathway,
                                          unsigned int model_len,
                                          unsigned int n_genes,
                                          bool has_passenger,
                                          double fbp,
                                          double bgp);
}

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
                           double *log_weights);

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
                        vector<shared_ptr<LinearProgressionParameters>> &ret_params);
#endif /* lpm_h */
