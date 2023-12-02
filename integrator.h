#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include<iostream>
#include"primitive_guassian.h"
#include<armadillo>
#include<cstdio>
#include<cmath>

//this file contains the integration functions for the trapezoid rule and the analytical overlap

namespace integrator {

//function to calculate integral by trapezoid rule
double trapezoid_rule(double (*func)(double), double lower_bound, double upper_bound, int num_intervals);

//function to calculate integral of the overlap of two primitive guassians by analytic solution 
double analytical_overlap_integral(int dimension, int l_a, int l_b, double alpha, double beta, 
                                        arma::vec R_a, arma::vec R_b, arma::vec R_p);
}

#endif