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

#include <spf/particle_genealogy.hpp>

#include "lpm_state.hpp"
#include "lpm_params.hpp"
#include "lpm_pg_proposal.hpp"

using namespace std;

bool path_exists(string path);
void compute_row_sum(const gsl_matrix &obs, vector<size_t> &row_sum);

// input related functions
gsl_matrix *get_first_n_lines(gsl_matrix *obs_matrix, size_t n);
gsl_matrix *read_data(string file_name, bool header = false);
gsl_matrix *read_csv(string file_name, bool header = false);
vector<size_t> *read_stages(string file_name);
vector<size_t> *read_ground_truth(string path_name);
unsigned int read_ground_truth_pathway_from_matrix(string file_name, unsigned int *ret); // returns number of pathways
void read_error_params(string path_name, double &fbp, double &bgp);
size_t read_true_model_length(string path_name);


// output related functions
void write_model_selection_output_to_file(string path, size_t n_patients,
                                          const vector<double> &log_evidences,
                                          const vector<double> &fbps,
                                          const vector<double> &bgps,
                                          const vector<vector<double>> &log_marginals);
void write_pg_output(string file,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<shared_ptr<LinearProgressionParameters> > &params);
void write_pg_output(string path,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<double> &log_marginal,
                     vector<shared_ptr<LinearProgressionParameters> > &params,
                     LPMParamProposal &pg_proposal);

void write_csv(string path, vector<unsigned int> data);
void write_matrix(string path, const gsl_matrix &data);

// data simulation functions
// sample a pathway into ret_pathway
// note: ret_pathway.size() == n_genes
void sample_pathway_uniform(gsl_rng *random, unsigned int n_pathways, vector<unsigned int> &ret_pathway);
// sample a pathway into ret_stages
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
