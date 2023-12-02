#include"integrator.h"
#include"util.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>

//function to integrate functions of one variable using the trapezoid rule
//inputs:
//      func() - any single variable function that returns a double
//      double lower_bound - the lower bound of the definite integral
//      double upper_bound - the upper bound of integration
//      int num_intervals - the number of sections to use in trapezoid rule
//returns:
//      double - the calculated value of the definite integral between the bounds
double integrator::trapezoid_rule(double (*func)(double), double lower_bound, double upper_bound, int num_intervals)
{
    double sum = 0;
    double h = (upper_bound - lower_bound)/num_intervals;

    for (int k = 1; k < num_intervals; k++)
    {
        sum += func(lower_bound + k*h);
    }

    sum = ((upper_bound - lower_bound)*(lower_bound + 2*sum) + func(lower_bound + num_intervals*h))/(2*num_intervals);

    return sum;

}

//function to calculate the analytical overlap integral in closed form
//inputs:
//      int dimension - the dimension (0, 1, 2) to integrate
//      int l_a - the value of l, m, or n of orbital A depending on the dimension chosen
//      int l_b - the value of l, m, or n of orbital B depending on the dimension chosen
//      double alpha - the exponent value of orbital A
//      double beta - the exponent value of orbital B
//      arma::vec R_b, R_b - the cartesian coordinates of the center of the orbitals A and B
//      arma::vec R_p - the center of the two orbitals as calculated by util::compute_center
//returns:
//      double - the calculated integral of the overlap between A and B from closed form solution
//
double integrator::analytical_overlap_integral(int dimension, int l_a, int l_b, double alpha, double beta, 
                                            arma::vec R_a, arma::vec R_b, arma::vec R_p)
{
    double X_a = R_a(dimension);
    double X_b = R_b(dimension);
    double X_p = R_p(dimension);

    double prefactor_term = std::exp(-alpha*beta*std::pow(X_a - X_b, 2) / (alpha + beta)) ;
    double square_root_term = std::pow(M_PI/(alpha+beta), 0.5);

    double summation_term = 0;

    for (int i = 0; i <=l_a; i++)
    {
        for (int j = 0; j <= l_b; j++)
        {
            if (((i + j) % 2) == 0)
            {
                summation_term += util::combination(l_a, i)*util::combination(l_b, j)
                * (util::double_factorial(i + j - 1) * std::pow(X_p - X_a, l_a - i) * std::pow(X_p - X_b, l_b - j)
                * (1.0 / std::pow(2*(alpha+beta), (i+j) / 2.0)));
            }
        }
    }

    double result = prefactor_term * square_root_term * summation_term;
    return result;
}

