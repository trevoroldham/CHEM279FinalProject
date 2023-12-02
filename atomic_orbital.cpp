#include"primitive_guassian.h"
#include"atomic_orbital.h"
#include"integrator.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>

//class containing the atomic_orbital object, defined as a guassian with l, m, n, alpha, and the center coordinates
atomic_orbital::atomic_orbital(arma::vec coords, int l, int m, int n, arma::vec alpha, arma::vec contraction_coefficient)
{

this->l = l;
this->n = n;
this->m = m;
this->coords = coords;
this->alpha = alpha;
this->contraction_coefficient = contraction_coefficient;
arma::vec norm(3);
this->normalization_constant = norm;


//initialize the three primitive guassians
primitive_guassian guassian_1 = primitive_guassian(coords, l, m, n, alpha(0), contraction_coefficient(0));
primitive_guassian guassian_2 = primitive_guassian(coords, l, m, n, alpha(1), contraction_coefficient(1));
primitive_guassian guassian_3 = primitive_guassian(coords, l, m, n, alpha(2), contraction_coefficient(2));

//save primitive guassians in a vector
guassians.push_back(guassian_1);
guassians.push_back(guassian_2);
guassians.push_back(guassian_3);

}

//function to calculate the normalization constants of each guassian
//inputs: none
//outputs: none
void atomic_orbital::set_normalization_constants()
{
    
    for (int i = 0; i < 3; i++)
    {
        double overlap;
        double overlap_x = integrator::analytical_overlap_integral(0, l, l, alpha(i), alpha(i), coords, coords, coords);
        double overlap_y = integrator::analytical_overlap_integral(1, m, m, alpha(i), alpha(i), coords, coords, coords);
        double overlap_z = integrator::analytical_overlap_integral(2, n, n, alpha(i), alpha(i), coords, coords, coords);
        overlap = overlap_x * overlap_y * overlap_z;

        this->normalization_constant(i) = 1.0/(std::pow(overlap, 0.5));
    }
};

