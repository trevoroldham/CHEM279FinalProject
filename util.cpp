#include"util.h"
#include"molecule.h"
#include"integrator.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>

//factorial: function to compute the factorial from a given integer
//inputs:
//      int n - the integer with which to calculate factorial
//outputs:
//      int - the factorial of n
int util::factorial(int n)
{

    if (n == 0 || n == 1)
    {
        return 1;
    }

    return n * factorial(n-1);
};

//double_factorial: function to compute the double factorial
//inputs:
//      int n - the integer from which to calculate double factorial
//outputs
//      int - the double factorial of n
int util::double_factorial(int n)
{
    if (n < 0)
    {
        return 1;
    }
    if (n == 0 || n == 1)
    {
        return 1;
    }

    return n*double_factorial(n-2);
    

}

//combination: function to compute the combinations of two integers
//inputs:
//      int m 
//      int n
//outputs:
//      int - the resulting number of combinations from m and n
int util::combination(int m, int n)
{
    if (n >= m)
    {
        return 1;
    }
    else {
    double result = util::factorial(m) / (util::factorial(n) * util::factorial(m - n));
    return result;
    }
    
}

//compute center: function to find the center of the interaction of two guassian orbital functions
//inputs:
//      arma::vec R_a - the cartesian coordinates center of the first orbital function
//      arma::vec R_b = the cartesian coordinates center of the second orbital function
//outputs:
//      aram::vec - the resulting center of the interaction
arma::vec util::compute_center(arma::vec R_a, arma::vec R_b, double alpha, double beta)
{
    arma::vec R_p = (alpha*R_a + beta*R_b) / (alpha + beta);
    return R_p;
}

//function to calculate the overlap integral between two primitive guassian functions
//inputs:
//      primitive_guassian a, b - two primitive guassian objects
//outputs:
//      double overlap - the overlap integral of the two guassians
double util::calculate_overlap(primitive_guassian a, primitive_guassian b)
{
    double overlap;
    double overlap_x;
    double overlap_y;
    double overlap_z;
    arma::vec center = util::compute_center(a.coords, b.coords, a.alpha, b.alpha);

 
    overlap_x = integrator::analytical_overlap_integral(0, a.l, b.l, a.alpha, b.alpha, a.coords, b.coords, center);
    overlap_y = integrator::analytical_overlap_integral(1, a.m, b.m, a.alpha, b.alpha, a.coords, b.coords, center);
    overlap_z = integrator::analytical_overlap_integral(2, a.n, b.n, a.alpha, b.alpha, a.coords, b.coords, center);
    overlap = overlap_x * overlap_y * overlap_z;  
    return overlap;
}

//function to calculate the overlap between two atomic orbital objects
//inputs:
//      atomic_orbital a, b - the two atomic orbital objects
//outputs:
//      double overlap - the overlap integral between a and b. this is the (a, b) element in the overlap matrix
double util::calculate_overlap(atomic_orbital a, atomic_orbital b)
{
    double overlap = 0;

    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            overlap += (a.contraction_coefficient(k) * b.contraction_coefficient(l) 
                        * a.normalization_constant(k) * b.normalization_constant(l)
                        * calculate_overlap(a.guassians[k], b.guassians[l]));
            
        }
    }
    return overlap;
}

//function to calculate the derivative of the overlap integral with respect to a nuclear perturbation
//used to calculate one row entry (3 x 1) in the derivative of the overlap matrix
//inputs:
//      atomic_orbital a, b - two atomic orbital objects
//outputs:
//      arma::vec result - (1 x 3) vector of the derivative in x, y, and z
arma::vec util::calculate_d_s_mu_nu(atomic_orbital a, atomic_orbital b)
{
    double dx = 0;
    double dy = 0;
    double dz = 0;

   
    arma::vec center = {0, 0, 0};
    arma::vec result = {0, 0, 0};

    int dim = 0;
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            center = util::compute_center(a.coords, b.coords, a.alpha(k), b.alpha(l));
            dx += (a.contraction_coefficient(k) * b.contraction_coefficient(l) 
                        * a.normalization_constant(k) * b.normalization_constant(l)
                        * util::calculate_d_sAB(0, a.guassians[k], b.guassians[l])
                        * integrator::analytical_overlap_integral(1, a.guassians[k].m, b.guassians[l].m, a.alpha(k), b.alpha(l), a.coords, b.coords, center)
                        * integrator::analytical_overlap_integral(2, a.guassians[k].n, b.guassians[l].n, a.alpha(k), b.alpha(l), a.coords, b.coords, center));
            
        }
    }

    
    dim = 1;
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            center = util::compute_center(a.coords, b.coords, a.alpha(k), b.alpha(l));
            dy += (a.contraction_coefficient(k) * b.contraction_coefficient(l) 
                        * a.normalization_constant(k) * b.normalization_constant(l)
                        * util::calculate_d_sAB(1, a.guassians[k], b.guassians[l])
                        * integrator::analytical_overlap_integral(0, a.guassians[k].l, b.guassians[l].l, a.alpha(k), b.alpha(l), a.coords, b.coords, center)
                        * integrator::analytical_overlap_integral(2, a.guassians[k].n, b.guassians[l].n, a.alpha(k), b.alpha(l), a.coords, b.coords, center));
            
        }
    }

    dim = 2;
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            center = util::compute_center(a.coords, b.coords, a.alpha(k), b.alpha(l));
            dz += (a.contraction_coefficient(k) * b.contraction_coefficient(l) 
                        * a.normalization_constant(k) * b.normalization_constant(l)
                        * util::calculate_d_sAB(2, a.guassians[k], b.guassians[l])
                        * integrator::analytical_overlap_integral(0, a.guassians[k].l, b.guassians[l].l, a.alpha(k), b.alpha(l), a.coords, b.coords, center)
                        * integrator::analytical_overlap_integral(1, a.guassians[k].m, b.guassians[l].m, a.alpha(k), b.alpha(l), a.coords, b.coords, center));
            
        }
    }

    result(0) = dx;
    result(1) = dy;
    result(2) = dz;

    return result;
}

//function to build the derivative of the overlap matrix
//inputs:
//      molecule mol - the input molecule
//outputs:
//      arma::field<arma::vec> s_mu_nu - (3 x n_basis x n_basis) field containing the derivative of the overlap wrt nuclear perturbation
arma::field<arma::vec> util::calculate_s_mu_nu_matrix(molecule mol)
{
    arma::field<arma::vec> s_mu_nu(mol.n_basis, mol.n_basis);

    for (int mu = 0; mu < mol.n_basis; mu++)
    {
        for (int nu = 0; nu < mol.n_basis; nu++)
        {
            if (mol.atom_map[mu] != mol.atom_map[nu])
            {
                s_mu_nu(mu, nu) = util::calculate_d_s_mu_nu(mol.orbitals[mu], mol.orbitals[nu]);
            }

            else
            {
                s_mu_nu(mu, nu) = {0, 0, 0};
            }
        }
    }

    return s_mu_nu;
}

//function to calculate the gamma value between two atomic orbitals
//inputs:
//      atomic_orbital a, b - the two atomic orbital objects
//outputs:
//      double result - the gamma value which is the (a, b) entry of gamma matrix
double util::calculate_gamma(atomic_orbital a, atomic_orbital b)
{
    arma::vec d_a_coefficients = a.contraction_coefficient;
    arma::vec d_b_coefficients = b.contraction_coefficient;

    arma::vec alpha_a_coefficients = a.alpha;
    arma::vec alpha_b_coefficients = b.alpha;


    //calculate d prime
    for (int i = 0; i < 3; i++)
    {
        d_a_coefficients(i) = d_a_coefficients(i) * a.normalization_constant(i);
        d_b_coefficients(i) = d_b_coefficients(i) * b.normalization_constant(i);
    }
    
    double result = 0;

    for (int k_1 = 0; k_1 < 3; k_1++)
    {
        for (int k_2 = 0; k_2 < 3; k_2++)
        {
            double sigma_a = 1.0 / (alpha_a_coefficients(k_1) + alpha_a_coefficients(k_2));

            for (int l_1 = 0; l_1 < 3; l_1++)
            {
                for (int l_2 = 0; l_2 < 3; l_2++)
                {
                    double sigma_b = 1.0 / (alpha_b_coefficients(l_1) + alpha_b_coefficients(l_2));
                    result += d_a_coefficients(k_1) * d_a_coefficients(k_2)
                                * d_b_coefficients(l_1) * d_b_coefficients(l_2)
                                * util::calculate_2e_integral(a, b, sigma_a, sigma_b);

                }
            }
        }
    }

    return result;
}

//function to calculate the two electron integral between two orbitals
//inputs:
//      atomic_orbital a, b - the two atomic orbitals
//      double sigma_a, sigma_b - the sigma values of the two atomic orbitals
//outputs:
//      double result - the integral value in electron volts
double util::calculate_2e_integral(atomic_orbital a, atomic_orbital b, double sigma_a, double sigma_b)
{
    double result=0;

    double u = std::pow(M_PI * sigma_a, 1.5) * std::pow(M_PI * sigma_b, 1.5);
    double rd = arma::norm(a.coords - b.coords, 2);
    double v2 = 1.0 / (sigma_a + sigma_b);
    double srT = std::sqrt(v2) * rd;
    
    if (rd < 1e-6)
    {
        result = u * std::sqrt(2*v2) * std::sqrt(2.0 / M_PI);
    }

    else 
    {
        result = (u / rd) * erf(srT);
    }

    return result*27.211;
}

//function to calculate [0]^0 term used in the function util::calculate_d_gamma_ab()
//inputs:
//      atomic_orbital a, b - two atomic orbital objects
//      double sigma_a, sigma_b - the sigma values calculated for a and b
//outputs:
//      arma::vec d_0_0 - (1 x 3) vector with the (x, y, z) derivatives of [0]^0 term in electron volts
arma::vec util::calculate_d_0_0(atomic_orbital a, atomic_orbital b, double sigma_a, double sigma_b)
{

    double u_ab = std::pow(M_PI * sigma_a, 1.5) * std::pow(M_PI * sigma_b, 1.5);
    double rd = arma::norm(a.coords - b.coords, 2);
    double v2 = 1.0 / (sigma_a + sigma_b);
    double srT = std::sqrt(v2) * rd;
    double T = std::pow(srT, 2);

    arma::vec d_0_0 = {0, 0, 0};

    if (rd != 0.0)
    {
        d_0_0 = (u_ab) * (1.0/std::pow(rd, 2)) * ((2*std::sqrt(v2) *exp(-T)) / std::sqrt(M_PI) - erf(srT)/rd) * (a.coords - b.coords);

    }
    return  d_0_0 * 27.211;
}

//function to calculate the derivative of an element of the gamma matrix wrt nuclear perturbation
//inputs:
//      int dim - the dimension along which to calculate
//      primitive_guassian a, b - two primitive guassians with which to calculate
//outputs:
//      double dx - the derivative of the dim dimension of the gamma value of (a, b)
double util::calculate_d_sAB(int dim, primitive_guassian a, primitive_guassian b)
{
    double dx = 0;

    arma::vec center = util::compute_center(a.coords, b.coords, a.alpha, b.alpha);
    int L = a.l + a.m + a.n;
    int l_a;
    int l_b;

    if (dim == 0)
    {
        l_a = a.l;
        l_b = b.l;
    }

    else if (dim == 1)
    {
        l_a = a.m;
        l_b = b.m;
    }

    else if (dim == 2)
    {
        l_a = a.n;
        l_b = b.n;
    }

    if (l_a == 0)
    {
        dx = 2 * a.alpha * integrator::analytical_overlap_integral(dim, l_a+1, l_b, a.alpha, b.alpha, a.coords, b.coords, center);
    }

    if (l_a == 1)
    {
        dx = (-1.0) * integrator::analytical_overlap_integral(dim, l_a - 1, l_b, a.alpha, b.alpha, a.coords, b.coords, center)
                + 2 * a.alpha * integrator::analytical_overlap_integral(dim, l_a + 1, l_b, a.alpha, b.alpha, a.coords, b.coords, center);
    }

    return dx;
}


//function to calculate the derivative of (a, b) element of the gamma matrix
//inputs:
//      atomic_orbital a, b - the two atomic orbital objects
//outputs:
//      arma::vec result - (3 x 1) vector with the (x, y, z) derivatives - the (a, b) entry of the field d_gamma_ab
arma::vec util::calculate_d_gamma_ab(atomic_orbital a, atomic_orbital b)
{
    arma::vec d_a_coefficients = a.contraction_coefficient;
    arma::vec d_b_coefficients = b.contraction_coefficient;

    arma::vec alpha_a_coefficients = a.alpha;
    arma::vec alpha_b_coefficients = b.alpha;


    //calculate d prime
    for (int i = 0; i < 3; i++)
    {
        d_a_coefficients(i) = d_a_coefficients(i) * a.normalization_constant(i);
        d_b_coefficients(i) = d_b_coefficients(i) * b.normalization_constant(i);
    }
    
    arma::vec result(3);

    for (int k_1 = 0; k_1 < 3; k_1++)
    {
        for (int k_2 = 0; k_2 < 3; k_2++)
        {
            double sigma_a = 1.0 / (alpha_a_coefficients(k_1) + alpha_a_coefficients(k_2));

            for (int l_1 = 0; l_1 < 3; l_1++)
            {
                for (int l_2 = 0; l_2 < 3; l_2++)
                {
                    double sigma_b = 1.0 / (alpha_b_coefficients(l_1) + alpha_b_coefficients(l_2));
                    result += d_a_coefficients(k_1) * d_a_coefficients(k_2)
                                * d_b_coefficients(l_1) * d_b_coefficients(l_2)
                                * util::calculate_d_0_0(a, b, sigma_a, sigma_b);

                }
            }
        }
    }

    return result;
}

//function to build the d_gamma_ab field (3 x n_atoms x n_atoms)
//inputs:
//      molecule mol - the given molecule
//outputs:
//      arma::field<arma::vec> result - the (3 x n_atoms x n_atoms) field containing the (x, y, z) derivative of the (a, b) element
arma::field<arma::vec> util::calculate_d_gamma_ab_matrix(molecule mol)
{
    arma::field<arma::vec> d_gamma_ab(mol.n_atoms, mol.n_atoms);

    int atom_a = 0;
    int orbital_a_index = 0;
    int atom_b = 0;
    int orbital_b_index = 0;
    arma::vec result;

    while (atom_a < mol.n_atoms)
    {
        atom_b  = 0;
        orbital_b_index = 0;

        while (atom_b < mol.n_atoms)
        {
            atomic_orbital a = mol.orbitals[orbital_a_index];
            atomic_orbital b = mol.orbitals[orbital_b_index];

            result = util::calculate_d_gamma_ab(a, b);
            d_gamma_ab(atom_a, atom_b) = result;
            
            if (mol.z_vals(atom_b) != 1)
            {
                atom_b++;
                orbital_b_index += 4;
            }

            else 
            {
                atom_b++;
                orbital_b_index++;
            }
           
        }

        if (mol.z_vals(atom_a) != 1)
        {
            atom_a++;
            orbital_a_index += 4;
        }

        else 
        {
            atom_a++;
            orbital_a_index++;
        }

    }

    return d_gamma_ab;

}


