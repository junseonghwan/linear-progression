//
//  data_util.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-02-01.
//

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "data_util.hpp"

// returns double array of length n_patients*n_gene
gsl_matrix *read_data(string file_name)
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
    getline (dat_file, line); // first line is just the header
    getline (dat_file, line);
    vector<string> results;
    boost::split(results, line, boost::is_any_of(","));
    fbp = stod(results[0]);
    bgp = stod(results[1]);
    dat_file.close();
}

size_t read_true_model_length(string path_name)
{
    return 0;
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

void write_model_selection_output_to_file(string file, const vector<double> &mean, const vector<double> &sd)
{
    ofstream myfile;
    myfile.open(file, ios::out);
    if (myfile.is_open()) {
        if (mean.size() != sd.size()) {
            cerr << "Error: vectors mean and sd are of different lengths. Mean vector size: " << mean.size() << ", SD vector size: " << sd.size() << endl;
            exit(-1);
        }
        // write data to file
        for (size_t l = 0; l < mean.size(); l++) {
            myfile << (l+1) << ", " << mean[l] << ", " << sd[l] << std::endl;
        }
        myfile.close();
    } else {
        cerr << "Error: cannot open " << file << endl;
        exit(-1);
    }
}
void write_pg_output(string path,
                     vector<ParticleGenealogy<LinearProgressionState> *> &states,
                     vector<LinearProgressionParameters *> &params)
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
