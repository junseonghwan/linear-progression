//
//  data_util.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-01.
//

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

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

void compute_row_sum(gsl_matrix &obs, vector<size_t> &row_sum)
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

