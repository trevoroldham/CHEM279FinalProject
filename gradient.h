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

namespace gradient {

arma::mat build_x_matrix(molecule mol);
arma::mat build_y_matrix(molecule mol);
arma::mat calculate_dv_dr(molecule mol);

void print_field(arma::field<arma::vec> field);

arma::mat calculate_gradient(molecule mol);
    
}

#endif