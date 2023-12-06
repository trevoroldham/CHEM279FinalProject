#ifndef SCF_H
#define SCF_H
#include"primitive_guassian.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>
#include"molecule.h"

namespace SCF
{
    void fixed_point_iteration(molecule & mol, double tolerance, bool verbose);
    void DIIS_FP(molecule & mol, double tolerance, int history, bool verbose);
    void DIIS(molecule & mol, double tolerance, int max_history, bool verbose);
}

#endif