#ifndef ATOMIC_ORBITAL_H
#define ATOMIC_ORBITAL_H
#include"primitive_guassian.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>

//class containing the atomic_orbital object, defined as a guassian with l, m, n, alpha, and the center coordinates
class atomic_orbital {
    public:
        int l;
        int m;
        int n;
        arma::vec alpha;
        arma::vec contraction_coefficient;
        arma::vec coords;
        arma::vec normalization_constant;
        

        //vector to store the three primitive guassian objects
        std::vector<primitive_guassian> guassians;

        //default constructor
        atomic_orbital(arma::vec coords, int l, int m, int n, arma::vec alpha, arma::vec contraction_coefficient);

        //function to set normalization constants
        void set_normalization_constants();
        
};

#endif