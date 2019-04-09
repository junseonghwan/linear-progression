//
//  data_util.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-01.
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <gsl/gsl_randist.h>

#include "data_util.hpp"

gsl_matrix *get_first_n_lines(gsl_matrix *obs_matrix, size_t n)
{
    gsl_matrix *ret = gsl_matrix_alloc(n, obs_matrix->size2);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < obs_matrix->size2; j++) {
            gsl_matrix_set(ret, i, j, gsl_matrix_get(obs_matrix, i, j));
        }
    }
    return ret;
}

bool path_exists(string path)
{
    return boost::filesystem::exists(path);
}

// returns double array of length n_patients*n_gene
gsl_matrix *read_data(string file_name, bool header)
{
    size_t n_patients = 0;
    size_t n_genes = 0;

    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }

    vector<string> results;
    vector<double> dat;
    
    if (header) {
        // skip the first line
        getline(dat_file, line);
    }
    
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(","));
        if (n_genes == 0) {
            n_genes = results.size();
        }
        else if (results.size() != n_genes) {
            cerr << "Error in the input file: " << file_name << endl;
            exit(-1);
        }
        for (size_t i = 0; i < n_genes; i++) {
            dat.push_back(stod(results[i]));
        }
        n_patients++;
    }
    dat_file.close();

    gsl_matrix *matrix = gsl_matrix_alloc(n_patients, n_genes);
    size_t idx = 0;
    for (size_t m = 0; m < n_patients; m++) {
        for (size_t n = 0; n < n_genes; n++) {
            gsl_matrix_set(matrix, m, n, dat[idx]);
            idx++;
        }
    }

    return matrix;
}

gsl_matrix *read_csv(string file_name, bool header)
{
    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }
    
    vector<string> results;
    vector<double> dat;
    
    if (header) {
        // skip the first line
        getline(dat_file, line);
    }
    
    size_t n_cols = 0;
    size_t n_rows = 0;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(","));
        if (n_cols == 0) {
            n_cols = results.size();
            if (n_cols == 0) {
                cerr << "Error: " << file_name << " contains 0 columns." << endl;
                exit(-1);
            }
        }
        else if (results.size() != n_cols) {
            cerr << "Error: Column numbers do not match in the input file: " << file_name << endl;
            exit(-1);
        }
        for (size_t i = 0; i < n_cols; i++) {
            dat.push_back(stod(results[i]));
        }
        n_rows++;
    }
    dat_file.close();
    
    gsl_matrix *matrix = gsl_matrix_alloc(n_rows, n_cols);
    size_t idx = 0;
    for (size_t m = 0; m < n_rows; m++) {
        for (size_t n = 0; n < n_cols; n++) {
            gsl_matrix_set(matrix, m, n, dat[idx]);
            idx++;
        }
    }
    
    return matrix;
}

vector<size_t> *read_stages(string file_name)
{
    string line;
    ifstream dat_file (file_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        return 0;
    }
    
    auto *stages = new vector<size_t>();
    while ( getline (dat_file, line) )
    {
        unsigned int stage = stoi(line);
        stages->push_back(stage);
    }
    dat_file.close();
    
    return stages;
}

vector<size_t> *read_ground_truth(string file_name)
{
    vector<size_t> *gt = new vector<size_t>();
    ifstream dat_file (file_name);
    string line;
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }
    vector<string> results;
    while ( getline (dat_file, line) )
    {
        boost::split(results, line, boost::is_any_of(","));
        //unsigned int gene = stoi(results[0]);
        unsigned int pathway = stoi(results[1]) - 1;
        gt->push_back(pathway);
    }
    dat_file.close();
    
    return gt;
}

unsigned int read_ground_truth_pathway_from_matrix(string file_name, unsigned int *ret)
{
    ifstream dat_file (file_name);
    string line;
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << file_name << endl;
        exit(-1);
    }

    size_t gene_idx;
    size_t pathway_idx;

    gsl_matrix *pathway_matrix = read_csv(file_name, false);
    for (gene_idx = 0; gene_idx < pathway_matrix->size1; gene_idx++) {
        for (pathway_idx = 0; pathway_idx < pathway_matrix->size1; pathway_idx++) {
            if (gsl_matrix_get(pathway_matrix, gene_idx, pathway_idx) == 1) {
                ret[gene_idx] = pathway_idx;
                break;
            }
        }
    }
    return pathway_matrix->size2;
}

void read_error_params(string path_name, double &fbp, double &bgp)
{
    string line;
    ifstream dat_file (path_name);
    if (!dat_file.is_open())
    {
        cerr << "Could not open the file: " << path_name << endl;
    }
    
    getline (dat_file, line); // first line is FBP
    fbp = stod(line);

    getline (dat_file, line); // second line is BGP
    bgp = stod(line);

    dat_file.close();

    
//    getline (dat_file, line); // first line is just the header
//    getline (dat_file, line);
//    vector<string> results;
//    boost::split(results, line, boost::is_any_of(","));
//    fbp = stod(results[0]);
//    bgp = stod(results[1]);
//    dat_file.close();
}

double *convert_to_array(vector<size_t> &pathway_membership, size_t n_pathways)
{
    size_t n_genes = pathway_membership.size();
    double *x = new double[n_genes * n_pathways];
    for (size_t n = 0; n < n_genes; n++) {
        size_t idx = pathway_membership[n];
        for (size_t k = 0; k < n_pathways; k++) {
            x[n * n_pathways + k] = idx == k ? 1 : 0;
        }
    }
    return x;
}

void write_model_selection_output_to_file(string path, size_t n_patients,
                                          const vector<double> &log_evidences,
                                          const vector<double> &fbps, const vector<double> &bgps,
                                          const vector<vector<double>> &log_marginals)
{
    size_t num_samples = fbps.size();
    size_t L = log_evidences.size();
    
    ofstream myfile;
    string evidence_file = path + "/model_evidences_" + to_string(n_patients) + ".csv";
    myfile.open(evidence_file, ios::out);
    if (myfile.is_open()) {
        // write log_evidences to a file
        for (size_t l = 0; l < L; l++) {
            myfile << (l+2) << ", " << log_evidences[l] << std::endl;
        }
        myfile.close();
    } else {
        cerr << "Error: cannot open " << evidence_file << endl;
        exit(-1);
    }
    
    string log_marginals_file = path + "/log_marginals_" + to_string(n_patients) + ".csv";
    myfile.open(log_marginals_file, ios::out);
    if (myfile.is_open()) {
        // write log_evidences to a file
        for (size_t i = 0; i < num_samples; i++) {
            for (size_t l = 0; l < L; l++) {
                myfile << (l+2) << ", " << fbps[i] << ", " << bgps[i] << ", " << log_marginals[l][i] << endl;
            }
        }
        myfile.close();
    } else {
        cerr << "Error: cannot open " << log_marginals_file << endl;
        exit(-1);
    }
    
}

void write_pg_output(string path,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<shared_ptr<LinearProgressionParameters> > &params)
{
    if (states.size() != params.size()) {
        cerr << "Error: vectors states and params are of different lengths. State vector size: " << states.size() << ", parameter vector size: " << params.size() << endl;
        exit(-1);
    }
    if (states.size() == 0 || params.size() == 0) {
        cout << "Warning: no posterior samples to write!" << endl;
        return;
    }
    
    string state_file_path = path + "/states.csv";
    string params_file_path = path + "/params.csv";

    ofstream state_file;
    ofstream params_file;
    state_file.open(state_file_path, ios::out);
    params_file.open(params_file_path, ios::out);
    if (!state_file.is_open()) {
        cerr << "Error: cannot open " << state_file_path << endl;
        exit(-1);
    }
    if (!params_file.is_open()) {
        cerr << "Error: cannot open " << params_file_path << endl;
        exit(-1);
    }

    // write data to file
    size_t len = states[0]->size();
    params_file << "FBP, BGP" << endl;
    for (size_t i = 0; i < states.size(); i++) {
        state_file << states[i]->get_state_at(len-1).to_string() << endl;
        params_file << params.at(i)->get_fbp() << ", " << params.at(i)->get_bgp() << endl;
    }

    state_file.close();
    params_file.close();
}

void write_pg_output(string path,
                     vector<shared_ptr<ParticleGenealogy<LinearProgressionState> > > &states,
                     vector<double> &log_marginal,
                     vector<shared_ptr<LinearProgressionParameters> > &params,
                     LPMParamProposal &pg_proposal)
{
    if (states.size() == 0) {
        cout << "Warning: no posterior samples to write!" << endl;
        return;
    }
    if (pg_proposal.get_bgps().size() != pg_proposal.get_fbps().size()) {
        cerr << "Error: number of BGPs sampled does not match the number of FBPs sampled." << endl;
        exit(-1);
    }

    string state_file_path = path + "/states.csv";
    string params_file_path = path + "/params.csv";
    string log_marginal_file_path = path + "/log_marginal.csv";

    ofstream state_file;
    ofstream params_file;
    ofstream log_marginal_file;
    state_file.open(state_file_path, ios::out);
    params_file.open(params_file_path, ios::out);
    log_marginal_file.open(log_marginal_file_path, ios::out);
    if (!state_file.is_open()) {
        cerr << "Error: cannot open " << state_file_path << endl;
        exit(-1);
    }
    if (!params_file.is_open()) {
        cerr << "Error: cannot open " << params_file_path << endl;
        exit(-1);
    }
    if (!log_marginal_file.is_open()) {
        cerr << "Error: cannot open " << log_marginal_file_path << endl;
        exit(-1);
    }

    // write states
    size_t len = states[0]->size();
    for (size_t i = 0; i < states.size(); i++) {
        state_file << states[i]->get_state_at(len-1).to_string() << endl;
    }
    state_file.close();

    // write all parameters
    const vector<double> &bgps = pg_proposal.get_bgps();
    const vector<double> &fbps = pg_proposal.get_fbps();
    params_file << "FBP, BGP" << endl;
    for (size_t i = 0; i < pg_proposal.get_bgps().size(); i++) {
        params_file << bgps[i] << ", " << fbps[i] << endl;
    }
    params_file.close();
    
    // write log marginals
    log_marginal_file << "FBP, BGP, LogMarginal" << endl;
    for (size_t i = 0; i < log_marginal.size(); i++) {
        log_marginal_file << params[i].get()->get_fbp() << ", " << params[i].get()->get_bgp() << ", " << log_marginal[i] << endl;
    }
    log_marginal_file.close();
}

void write_csv(string path, vector<unsigned int> data)
{
    ofstream f;
    f.open(path, ios::out);
    for (size_t i = 0; i < data.size(); i++) {
        f << data[i] << endl;
    }
    f.close();
}

void write_matrix(string path, const gsl_matrix &data)
{
    ofstream f;
    f.open(path, ios::out);
    if (f.is_open()) {
        for (size_t i = 0; i < data.size1; i++) {
            for (size_t j = 0; j < data.size2; j++) {
                f << gsl_matrix_get(&data, i, j);
                if (j < data.size2 - 1) {
                    f << ", ";
                }
            }
            f << "\n";
        }
        f.close();
    } else {
        cerr << "Error: cannot open file: " << path << endl;
        exit(-1);
    }
}

void compute_row_sum(const gsl_matrix &obs, vector<size_t> &row_sum)
{
    // compute the row sum for obs matrix
    size_t M = obs.size1;
    size_t N = obs.size2;
    if (row_sum.size() != M) {
        cerr << "Error: dimension for row_sum: " << row_sum.size() << " does not num rows of obs: " << M << endl;
        exit(-1);
    }
    for (size_t m = 0; m < M; m++) {
        size_t sum_m = 0;
        for (size_t n = 0; n < N; n++) {
            sum_m += (size_t)gsl_matrix_get(&obs, m, n);
        }
        row_sum[m] = sum_m;
    }
}

void sample_pathway_uniform(gsl_rng *random, unsigned int num_pathways, vector<unsigned int> &ret_pathway)
{
    size_t n_genes = ret_pathway.size();
    size_t num_divider_locations = n_genes - 1;
    size_t num_dividers = num_pathways;

    size_t *possible_divider_positions = new size_t[num_divider_locations];
    size_t *dividers = new size_t[num_pathways];
    size_t *genes = new size_t[n_genes];

    // populate divider locations
    for (size_t n = 0; n < num_divider_locations; n++) {
        possible_divider_positions[n] = n + 1;
    }

    // populate genes
    for (size_t n = 0; n < n_genes; n++) {
        genes[n] = n;
    }
    
    // select divider locations
    gsl_ran_choose(random, dividers, num_dividers - 1, possible_divider_positions, num_divider_locations, sizeof(size_t));
    // last divider location is the last index
    dividers[num_dividers - 1] = num_divider_locations + 1;
    sort(dividers, dividers + num_dividers);

    gsl_ran_shuffle(random, genes, n_genes, sizeof(size_t));
    size_t n = 0;
    for (size_t k = 0; k < num_pathways; k++) {
        for (; n < n_genes; n++) {
            size_t g = genes[n];
            if (n < dividers[k]) {
                ret_pathway[g] = k;
            } else {
                break;
            }
        }
    }
    
    delete [] genes;
    delete [] possible_divider_positions;
    delete [] dividers;
}

void sample_stages_uniform(gsl_rng *random, unsigned int n_pathways, vector<unsigned int> &ret_stages)
{
    unsigned int n_patients = ret_stages.size();
    unsigned int i;
    unsigned int stage;
    for (i = 0; i < n_patients; i++) {
        stage = gsl_rng_uniform_int(random, n_pathways);
        ret_stages[i] = stage;
    }
}

gsl_matrix *simulate_data(gsl_rng *random,
                          unsigned int n_pathways,
                          const vector<unsigned int> &pathway,
                          const vector<unsigned int> &patient_stages)
{
    unsigned int n_patients = patient_stages.size();
    unsigned int n_genes = pathway.size();
    unsigned int k, g, n;

    gsl_matrix *data_matrix = gsl_matrix_alloc(n_patients, n_genes);

    // convert pathway to map : k -> vector
    unordered_map<unsigned int, vector<unsigned int>> pathway2genes;
    for (g = 0; g < n_genes; g++) {
        k = pathway[g];
        pathway2genes[k].push_back(g);
    }
    
    for (n = 0; n < n_patients; n++) {
        for (g = 0; g < n_genes; g++) {
            gsl_matrix_set(data_matrix, n, g, 0);
        }
    }

    for (n = 0; n < n_patients; n++) {
        for (k = 0; k <= patient_stages[n]; k++) {
            // sample gene from pathway2genes[k]
            g = pathway2genes[k][gsl_rng_uniform_int(random, pathway2genes[k].size())];
            gsl_matrix_set(data_matrix, n, g, 1.0);
        }
    }
    return data_matrix;
}

void add_noise(gsl_rng *random,
               double fbp,
               double bgp,
               gsl_matrix *matrix)
{
    size_t i, j;
    double val, unif;
    for (i = 0; i < matrix->size1; i++) {
        for (j = 0; j < matrix->size2; j++) {
            val = gsl_matrix_get(matrix, i, j);
            unif = gsl_ran_flat(random, 0.0, 1.0);
            if (val == 0) {
                // set val to 1 (from 0) with prob bgp
                if (unif < bgp) {
                    //cout << "Background mutation at (" << i << ", " << j << ")" << endl;
                    val = 1;
                }
            } else {
                // set val to 0 (from 1) with prob fbp
                if (unif < fbp) {
                    //cout << "Flip back at (" << i << ", " << j << ")" << endl;
                    val = 0;
                }
            }
            gsl_matrix_set(matrix, i, j, val);
        }
    }
}
