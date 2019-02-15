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

using namespace std;

gsl_matrix *read_data(string file_name, bool header = false);
vector<size_t> *read_stages(string file_name);
vector<size_t> *read_ground_truth(string path_name);
void read_error_params(string path_name, double &fbp, double &bgp);
size_t read_true_model_length(string path_name);
void write_model_selection_output_to_file(string file, const vector<double> &mean, const vector<double> &sd);
void write_pg_output(string file,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<shared_ptr<LinearProgressionParameters> > &params);

#endif /* data_util_hpp */
