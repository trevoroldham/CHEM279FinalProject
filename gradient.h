#ifndef GRADIENT_H
#define GRADIENT_H
#include"primitive_guassian.h"
#include"atomic_orbital.h"
#include"molecule.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>

//namespace to hold information and functions for the calculation of the CNDO/2 analytical gradient

namespace gradient {

//functions to calculate relevant terms in the analytical gradient expression
arma::mat build_x_matrix(molecule mol);
arma::mat build_y_matrix(molecule mol);
arma::mat calculate_dv_dr(molecule mol);

//function to print a field object to std::cout
void print_field(arma::field<arma::vec> field);

//function to calculate gradient from CNDO/2 model
arma::mat calculate_gradient(molecule mol);
    
}

#endif