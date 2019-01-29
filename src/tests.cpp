//
//  tests.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-28.
//

#include "tests.hpp"

//bool test_likelihood1()
//{
//    double fbp = 0.072;
//    double bgp = 0.0321;
//
//    int n_patients = 2;
//    int n_pathways = 3;
//    int n_genes = 3;
//
//    double y[] = {
//        1, 0, 1,
//        1, 1, 1,
//    };
//
//    double x[] = {
//        1, 0, 0,
//        0, 1, 0,
//        0, 0, 1,
//    };
//
//    double r[n_patients*n_pathways];
//    int pathway_sizes[] = {1, 1, 1};
//
//    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
//    gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
//    gsl_matrix_view R = gsl_matrix_view_array(r, n_patients, n_pathways);
//
//    /* Compute C = A B */
//    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
//                    1.0, &Y.matrix, &X.matrix,
//                    0.0, &R.matrix);
//
//
//    // test likelihood calculation
//    double truth[] = {-3.546249, -0.2241706};
//
//    int stages[] = {0, 2};
//    for (int m = 0; m < n_patients; m++) {
//        double log_lik = 0.0;
//        for (int k = 0; k < n_pathways; k++) {
//            double val = gsl_matrix_get(&R.matrix, m, k);
//            if (k <= stages[m]) {
//                log_lik += compute_log_lik(val, pathway_sizes[k], bgp, fbp);
//            } else {
//                log_lik += compute_log_lik_bg(val, pathway_sizes[k], bgp, fbp);
//            }
//        }
//        if (abs(truth[m] - log_lik) > TOL) {
//            cerr << "Likelihood test 1: error in likelihood computation for sample " << m << endl;
//            cerr << "Expected: " << truth[m] << endl;
//            cerr << "Obtained: " << log_lik << endl;
//            return false;
//        }
//    }
//    cout << "Likelihood test 1 passed!" << endl;
//    return true;
//}
//
//bool test_likelihood2()
//{
//    // simple test:
//    // read in the data
//    // read in the true pathway
//    // compute the likelihood
//
//    double fbp = 0.027;
//    double bgp = 0.1321;
//
//    int n_patients = 3;
//    int n_pathways = 3;
//    int n_genes = 7;
//
//    double y[] = {
//        1, 0, 1, 0, 0, 0, 0,
//        1, 1, 1, 0, 1, 1, 1,
//        0, 0, 1, 1, 0, 0, 0
//    };
//
//    double x[] = {
//        1, 0, 0,
//        1, 0, 0,
//        0, 1, 0,
//        0, 1, 0,
//        0, 1, 0,
//        0, 1, 0,
//        0, 0, 1
//    };
//
//    double r[n_patients*n_pathways];
//    int pathway_sizes[] = {2, 4, 1};
//
//    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
//    gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
//    gsl_matrix_view R = gsl_matrix_view_array(r, n_patients, n_pathways);
//
//    /* Compute C = A B */
//    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
//                    1.0, &Y.matrix, &X.matrix,
//                    0.0, &R.matrix);
//
//
//    // test likelihood calculation
//    double truth[] = {-4.21657, -4.503214, -4.839183};
//
//    int stages[] = {2, 2, 1};
//    for (int m = 0; m < n_patients; m++) {
//        double log_lik = 0.0;
//        for (int k = 0; k < n_pathways; k++) {
//            double val = gsl_matrix_get(&R.matrix, m, k);
//            if (k <= stages[m]) {
//                log_lik += compute_log_lik(val, pathway_sizes[k], bgp, fbp);
//            } else {
//                log_lik += compute_log_lik_bg(val, pathway_sizes[k], bgp, fbp);
//            }
//        }
//        if (abs(truth[m] - log_lik) > TOL) {
//            cerr << "Likelihood test 2: error in likelihood computation for sample " << m << endl;
//            cerr << "Expected: " << truth[m] << endl;
//            cerr << "Obtained: " << log_lik << endl;
//            return false;
//        }
//    }
//    cout << "Likelihood test 2 passed!" << endl;
//    return true;
//}
//
//bool test_likelihood3()
//{
//    // simple test:
//    // read in the data
//    // read in the true pathway
//    // compute the likelihood
//
//    double fbp = 0.027;
//    double bgp = 0.1321;
//
//    int n_patients = 1;
//    int n_pathways = 2;
//    int n_genes = 4;
//
//    double y[] = {
//        1, 0, 1, 0
//    };
//
//    vector<double> pathway;
//    pathway.push_back(0);
//    pathway.push_back(0);
//    pathway.push_back(1);
//    pathway.push_back(1);
//
//    double x[] = {
//        1, 0,
//        1, 0,
//        0, 1,
//        0, 1
//    };
//
//    double r[n_patients*n_pathways];
//
//    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
//    gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
//    gsl_matrix_view R = gsl_matrix_view_array(r, n_patients, n_pathways);
//
//    /* Compute C = A B */
//    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
//                    1.0, &Y.matrix, &X.matrix,
//                    0.0, &R.matrix);
//
//    LinearProgressionState true_state(n_genes, n_pathways);
//    for (size_t g = 0; g < n_genes; g++) {
//        true_state.update_pathway_membership(g, pathway[g]);
//    }
//
//    LinearProgressionParameters params(fbp, bgp);
//
//    // test likelihood calculation where patient stage is marginalized
//    double truth = -0.2028664;
//    double log_lik = compute_pathway_likelihood(Y, true_state, params);
//    if (abs(truth - log_lik) > TOL) {
//        cerr << "Likelihood test 3: error in the likelihood computation marginalization over stages." << endl;
//        cerr << "Expected: " << truth << endl;
//        cerr << "Obtained: " << log_lik << endl;
//        return false;
//    }
//    return true;
//}

//bool test_lpm_model_importance_sampling()
//{
//    // for small number of genes and pathway, we can compute the exact value of the marginal likelihood by
//    // enumerating over all possibilities
//    // then, we can check that the marginal likelihood estimate from SMC approaches this value
//
//    double fbp = 0.027;
//    double bgp = 0.1321;
//
//    int n_patients = 3;
//    int n_pathways = 2;
//    int n_genes = 3;
//
//    double y[] = {
//        1, 1, 0,
//        1, 0, 1,
//        1, 1, 1
//    };
//    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
//    vector<size_t> stages(n_patients, 0);
//    stages[0] = 1;
//    stages[1] = 0;
//    stages[2] = 2;
//
//    // enumerate over possible values of x using 3 for loops
//    double exact_log_lik = DOUBLE_NEG_INF;
//    vector<size_t> pathway_membership(n_genes, 0);
//    for (size_t i = 0; i < n_pathways; i++) {
//        pathway_membership.at(0) = i;
//        for (size_t j = 0; j < n_pathways; j++) {
//            pathway_membership.at(1) = j;
//            for (size_t k = 0; k < n_pathways; k++) {
//                pathway_membership.at(2) = k;
//
//                // convert pathway membership to matrix x
//                double *x = convert_to_array(pathway_membership, n_pathways);
//                gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
//                double log_val = compute_pathway_likelihood(Y, X, stages, fbp, bgp);
//                exact_log_lik = log_add(exact_log_lik, log_val);
//                cout << "===" << endl;
//                cout << log_val << ", " << exact_log_lik << endl;
//                for (size_t l = 0; l < pathway_membership.size(); l++)
//                    cout << pathway_membership[l] << endl;
//                cout << "===" << endl;
//            }
//        }
//    }
//    exact_log_lik += log(1./8); // multiply by the prior
//    cout << "Exact marginal likelihood: " << exact_log_lik << endl;
//
//    // run importance sampling with proposal from the prior, to get an estimate of the marginal log likelihood
//    SMCOptions smc_options;
//    smc_options.ess_threshold = 1.0;
//    smc_options.num_particles = 100000;
//    smc_options.resample_last_round = false;
//
//    size_t n_smc_iter = 1;
//    bool allocate_passenger_pathway = false;
//    LinearProgressionModel *model = new LinearProgressionModel(n_genes, n_pathways, n_smc_iter, 10, &Y, allocate_passenger_pathway);
//    LinearProgressionParameters params(fbp, bgp);
//    params.set_patient_progression_stages(stages);
//
//    SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
//    smc.run_smc(params);
//    ParticlePopulation<LinearProgressionState> *pop = smc.get_curr_population();
//    double logZ = pop->get_log_norm();
//    cout << "Estimated logZ: " << logZ << endl;
//    if (abs(exp(exact_log_lik) - exp(logZ)) > TOL) {
//        cerr << "Importance sampling estimate is incorrect." << endl;
//        return false;
//    }
//    return true;
//}
//
//bool test_lpm_model_smc()
//{
//    double fbp = 0.027;
//    double bgp = 0.1321;
//
//    int n_patients = 3;
//    int n_pathways = 2;
//    int n_genes = 3;
//
//    double y[] = {
//        1, 1, 0,
//        1, 0, 1,
//        1, 1, 1
//    };
//    gsl_matrix_view Y = gsl_matrix_view_array(y, n_patients, n_genes);
//    vector<size_t> stages(n_patients, 0);
//    stages[0] = 1;
//    stages[1] = 0;
//    stages[2] = 2;
//
//    // enumerate over possible values of x using 3 for loops
//    double exact_log_lik = DOUBLE_NEG_INF;
//    vector<size_t> pathway_membership(n_genes, 0);
//    for (size_t i = 0; i < n_pathways; i++) {
//        pathway_membership.at(0) = i;
//        for (size_t j = 0; j < n_pathways; j++) {
//            pathway_membership.at(1) = j;
//            for (size_t k = 0; k < n_pathways; k++) {
//                pathway_membership.at(2) = k;
//
//                // convert pathway membership to matrix x
//                double *x = convert_to_array(pathway_membership, n_pathways);
//                gsl_matrix_view X = gsl_matrix_view_array(x, n_genes, n_pathways);
//                double log_val = compute_pathway_likelihood(Y, X, stages, fbp, bgp);
//                exact_log_lik = log_add(exact_log_lik, log_val);
//                cout << "===" << endl;
//                cout << log_val << ", " << exact_log_lik << endl;
//                for (size_t l = 0; l < pathway_membership.size(); l++)
//                    cout << pathway_membership[l] << endl;
//                cout << "===" << endl;
//            }
//        }
//    }
//    exact_log_lik += log(1./8);
//    cout << "Exact marginal likelihood: " << exact_log_lik << endl;
//
//    // run importance sampling with proposal from the prior, to get an estimate of the marginal log likelihood
//    SMCOptions smc_options;
//    smc_options.ess_threshold = 0.8;
//    smc_options.num_particles = 2000;
//    smc_options.resample_last_round = false;
//    smc_options.track_population = true;
//    size_t n_smc_iter = 200;
//
//    bool allocate_passenger_pathway = false;
//    LinearProgressionModel *model = new LinearProgressionModel(n_genes, n_pathways, n_smc_iter, 5, &Y, allocate_passenger_pathway);
//    LinearProgressionParameters params(fbp, bgp);
//    params.set_patient_progression_stages(stages);
//
//    SMC<LinearProgressionState, LinearProgressionParameters> smc(model, &smc_options);
//    smc.run_smc(params);
//    double logZ = smc.get_log_marginal_likelihood();
//    cout << "=====" << endl;
//    cout << "Estimated logZ: " << logZ << endl;
//    if (abs(exp(exact_log_lik) - exp(logZ)) > TOL) {
//        cerr << "SMC estimate is incorrect." << endl;
//        return false;
//    }
//    return true;
//}
//
//void run_tests()
//{
//    if (!test_likelihood1()) {
//        cerr << "Error: likelihood test 1" << endl;
//        exit(-1);
//    }
//    if (!test_likelihood2()) {
//        cerr << "Error: likelihood test 2" << endl;
//        exit(-1);
//    }
//    if (!test_lpm_model_importance_sampling()) {
//        cerr << "Error: importance sampling" << endl;
//        exit(-1);
//    }
//    if (!test_lpm_model_smc()) {
//        cerr << "Error: smc sampler" << endl;
//        exit(-1);
//    }
//}

