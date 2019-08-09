/*
 * numerical_utils.cpp
 *
 *  Created on: 7 Mar 2018
 *      Author: seonghwanjun
 */

#include "numerical_utils.hpp"

#include <cmath>
#include <gsl/gsl_sf.h>

const double DOUBLE_INF = std::numeric_limits<double>::infinity();
const double DOUBLE_NEG_INF = -std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

double log_binomial_pdf(const unsigned int k, const double p, const unsigned int n)
{
    if (k > n)
    {
        return 0;
    }
    else
    {
        double logP;
        
        if (p == 0)
        {
            logP = (k == 0) ? 0 : DOUBLE_NEG_INF;
        }
        else if (p == 1)
        {
            logP = (k == n) ? 0 : DOUBLE_NEG_INF;
        }
        else
        {
            double ln_Cnk = gsl_sf_lnchoose (n, k);
            logP = ln_Cnk + k * log (p) + (n - k) * log1p(-p);
        }

        return logP;
    }
}

double normalize(const vector<double> &log_weights, vector<double> &weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    normalize(log_weights, weights, log_norm);
    return log_norm;
}

void normalize(const vector<double> &log_weights, vector<double> &weights, double log_norm)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        weights[i] = exp(log_weights[i] - log_norm);
        sum += weights[i];
    }
    if (abs(sum - 1.0) > 1e-4) {
        cerr << "Error in normalization. Check that log_weights and log_norm are correctly calculated." << endl;
        cerr << ceil(sum*100000)/100000.0 << " != 1.0" << endl;
        cerr << "log_norm: " << log_norm << endl;
        exit(-1);
    }
}

void normalize_destructively(vector<double> &log_weights, double log_norm)
{
    // compute the lognorm
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        log_weights[i] = exp(log_weights[i] - log_norm);
    }
}

double normalize_destructively(vector<double> &log_weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    normalize_destructively(log_weights, log_norm);
    return log_norm;
}

double normalize_destructively(double *log_weights, int size)
{
	// compute the lognorm
	double log_norm = log_add(log_weights, size);
	for (int i = 0; i < size; i++)
	{
		log_weights[i] = exp(log_weights[i] - log_norm);
	}
    return log_norm;
}

double log_add(double x, double y)
{
	// make x the max
	if (y > x) {
		double temp = x;
		x = y;
		y = temp;
	}
	// now x is bigger
	if (x == DOUBLE_NEG_INF) {
		return x;
	}
	double negDiff = y - x;
	if (negDiff < -20) {
		return x;
	}
	return x + log(1.0 + exp(negDiff));
}

// is this useful? or even, make sense?
double log_subtract(double x, double y)
{
    if (x < y) { // (log(e^x - e^y) = log(neg number) = -Inf
        return DOUBLE_NEG_INF;
    }
    double negDiff = y - x;
    if (negDiff < -20) {
        return x;
    }
    return x + log(1.0 - exp(negDiff));
}


double log_add(vector<double> x)
{
    double max = DOUBLE_NEG_INF;
    double maxIndex = 0;
    for (unsigned int i = 0; i < x.size(); i++)
    {
        if (x[i] > max) {
            max = x[i];
            maxIndex = i;
        }
    }
    if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
    // compute the negative difference
    double threshold = max - 20;
    double sumNegativeDifferences = 0.0;
    for (unsigned int i = 0; i < x.size(); i++) {
        if (i != maxIndex && x[i] > threshold) {
            sumNegativeDifferences += exp(x[i] - max);
        }
    }
    if (sumNegativeDifferences > 0.0) {
        return max + log(1.0 + sumNegativeDifferences);
    } else {
        return max;
    }
    
}

double log_add(double *x, int size)
{
	double max = DOUBLE_NEG_INF;
	double maxIndex = 0;
	for (int i = 0; i < size; i++)
	{
		if (x[i] > max) {
			max = x[i];
			maxIndex = i;
		}
	}
	if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
	  // compute the negative difference
	  double threshold = max - 20;
	  double sumNegativeDifferences = 0.0;
	  for (int i = 0; i < size; i++) {
	    if (i != maxIndex && x[i] > threshold) {
	      sumNegativeDifferences += exp(x[i] - max);
	    }
	  }
	  if (sumNegativeDifferences > 0.0) {
	    return max + log(1.0 + sumNegativeDifferences);
	  } else {
	    return max;
	  }

}

void add(double *x, double c, size_t size)
{
    for (size_t s = 0; s < size; s++)
    {
        x[s] += c;
    }
}

void multiply(double *x, double c, double *ret, size_t size)
{
    for (size_t s = 0; s < size; s++)
    {
        ret[s] = x[s] * c;
    }
}

template <typename T>
void print_vector(const vector<T> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
}

