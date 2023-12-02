#include"primitive_guassian.h"
#include"atomic_orbital.h"
#include"molecule.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>

namespace util {

    int factorial(int n);
    int double_factorial(int n);
    int combination(int m, int n);

    arma::vec compute_center(arma::vec R_a, arma::vec R_b, double alpha, double beta);

    double calculate_overlap(primitive_guassian a, primitive_guassian b);
    double calculate_overlap(atomic_orbital a, atomic_orbital b);
    double calculate_gamma(atomic_orbital a, atomic_orbital b);
    double calculate_2e_integral(atomic_orbital a, atomic_orbital b, double sigma_a, double sigma_b);
    arma::vec calculate_d_0_0(atomic_orbital a, atomic_orbital b, double sigma_a, double sigma_b);
    arma::vec calculate_d_gamma_ab(atomic_orbital a, atomic_orbital);
    arma::field<arma::vec> calculate_d_gamma_ab_matrix(molecule mol);
    double calculate_d_sAB(int dim, primitive_guassian a, primitive_guassian b);
    arma::vec calculate_d_s_mu_nu(atomic_orbital a, atomic_orbital b);
    arma::field<arma::vec> calculate_s_mu_nu_matrix(molecule mol);

}