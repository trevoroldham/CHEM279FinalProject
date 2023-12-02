#ifndef PRIMITIVE_GUASSIAN_H
#define PRIMITIVE_GUASSIAN_H
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>

//class containing the atomic_orbital object, defined as a guassian with l, m, n, alpha, and the center coordinates
class primitive_guassian {
    public:
        int l;
        int m;
        int n;
        double alpha;
        double contraction_coefficient;
        arma::vec coords;

        primitive_guassian(arma::vec coords, int l, int m, int n, double alpha, double contraction_coefficient);

        //std::vec<guassian> g(arma::vec coords, int l, int m, int n, double alpha, double contraction_coefficient)
};

#endif