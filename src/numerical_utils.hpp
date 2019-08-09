/*
 * numerical_utils.hpp
 *
 *  Created on: 6 Mar 2018
 *      Author: seonghwanjun
 */

#ifndef SRC_NUMERICAL_UTILS_HPP_
#define SRC_NUMERICAL_UTILS_HPP_

#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

extern const double DOUBLE_INF;
extern const double DOUBLE_NEG_INF;
extern const double NaN;

// returns log sum as by-product
double log_binomial_pdf(const unsigned int k, const double p, const unsigned int n);
double normalize(const vector<double> &log_weights, vector<double> &weights);
void normalize(const vector<double> &log_weights, vector<double> &weights, double log_norm);
void normalize_destructively(vector<double> &log_weights, double log_norm);
double normalize_destructively(vector<double> &log_weights);
double normalize_destructively(double *log_weights, int size);
double log_add(double x, double y);
double log_subtract(double x, double y);
double log_add(double *x, int size);
double log_add(vector<double> x);

void add(double *x, double c, size_t size);
void multiply(double *x, double c, double *ret, size_t size);

template <typename T>
void print_vector(const vector<T> &v);

#endif /* SRC_NUMERICAL_UTILS_HPP_ */
