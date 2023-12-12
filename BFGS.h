#ifndef BFGS_H
#define BFGS_H
#include"primitive_guassian.h"
#include"atomic_orbital.h"
#include"molecule.h"
#include"gradient.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>

//namespace to hold information and functions for the calculation of the CNDO/2 analytical gradient

namespace BFGS 
{
    void minimize(molecule & mol, double tolerance, int max_iter);
    double line_search(molecule mol, arma::vec p, arma::vec grad, int max_iter);

}

#endif