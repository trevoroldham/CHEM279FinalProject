#include"primitive_guassian.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>

primitive_guassian::primitive_guassian(arma::vec coords, int l, int m, int n, double alpha, double contraction_coefficient)
{
    this->l = l;
    this->m = m;
    this->n = n;
    this->alpha = alpha;
    this->contraction_coefficient = contraction_coefficient;
    this->coords = coords;
}