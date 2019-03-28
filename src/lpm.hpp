//
//  lpm.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-03-27.
//

#ifndef lpm_h
#define lpm_h

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
}

#endif /* lpm_h */
