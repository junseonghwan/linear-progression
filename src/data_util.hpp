//
//  data_util.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-01.
//

#ifndef data_util_hpp
#define data_util_hpp

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>

#include <spf/pg.hpp>

#include "lpm_state.hpp"
#include "lpm_params.hpp"
#include "lpm_pg_proposal.hpp"

using namespace std;

bool path_exists(string path);
void compute_row_sum(const gsl_matrix &obs, vector<unsigned int> &row_sum);

// input related functions
gsl_matrix *get_first_n_lines(gsl_matrix *obs_matrix, unsigned int n);
gsl_matrix *read_data(string file_name, string sep, bool header = false);
gsl_matrix *read_csv(string file_name, bool header = false);
vector<unsigned int> *read_stages(string file_name);
vector<unsigned int> *read_ground_truth(string path_name);
unsigned int read_ground_truth_pathway_from_matrix(string file_name, unsigned int *ret); // returns number of pathways
void read_error_params(string path_name, double &fbp, double &bgp);
unsigned int read_true_model_length(string path_name);

// output related functions
void write_model_selection_output_to_file(string path, unsigned int n_patients,
                                          const vector<double> &log_evidences,
                                          const vector<double> &fbps,
                                          const vector<double> &bgps,
                                          const vector<vector<double>> &log_marginals);
void write_pg_output(string path,
                     ParticleGibbs<LinearProgressionState, LinearProgressionParameters> &pg,
                     LPMParamProposal &pg_proposal);
void write_vector(string path, vector<unsigned int> data);
void write_matrix_as_csv(string path, const gsl_matrix &data);

// data simulation functions
// sample a pathway into ret_pathway
// note: ret_pathway.size() == n_genes
void sample_pathway_from_prior(gsl_rng *random, unsigned int n_pathways, vector<unsigned int> &ret_pathway);
// sample a patient progression stage into ret_stages
// note: ret_stages.size() == n_patients
void sample_stages_uniform(gsl_rng *random, unsigned int n_pathways, vector<unsigned int> &ret_stages);
// generate the clean matrix (no noise)
gsl_matrix *simulate_data(gsl_rng *random,
                          unsigned int n_pathways,
                          const vector<unsigned int> &pathway,
                          const vector<unsigned int> &patient_stages);
// add noise to the matrix (directly)
void add_noise(gsl_rng *random,
               double fbp,
               double bgp,
               gsl_matrix *matrix);

#endif /* data_util_hpp */
