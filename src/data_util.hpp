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

gsl_matrix *get_first_n_lines(gsl_matrix *obs_matrix, size_t n);
bool path_exists(string path);
gsl_matrix *read_data(string file_name, bool header = false);
gsl_matrix *read_csv(string file_name, bool header);
vector<size_t> *read_stages(string file_name);
vector<size_t> *read_ground_truth(string path_name);
void read_error_params(string path_name, double &fbp, double &bgp);
size_t read_true_model_length(string path_name);
void write_model_selection_output_to_file(string path, size_t n_patients,
                                          const vector<double> &log_evidences,
                                          const vector<double> &fbps, const vector<double> &bgps,
                                          const vector<vector<double>> &log_marginals);
void write_pg_output(string file,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<shared_ptr<LinearProgressionParameters> > &params);
void write_pg_output(string path,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<double> &log_marginal,
                     vector<shared_ptr<LinearProgressionParameters> > &params,
                     LPMParamProposal &pg_proposal);
void compute_row_sum(gsl_matrix &obs, vector<size_t> &row_sum);

#endif /* data_util_hpp */
