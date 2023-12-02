#include"util.h"
#include"gradient.h"
#include"molecule.h"
#include"integrator.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>

//function to build the x matrix in the CNDO/2 analytical gradient expression
//inputs:
//      molecule mol - the molecule object with which to calculate the gradient
//outputs:
//      arma::mat - a matrix of dimension (n_basis x n_basis)
arma::mat gradient::build_x_matrix(molecule mol)
{
    arma::mat x(mol.n_basis, mol.n_basis);
    int z_a;
    int z_b;

    for (int mu = 0; mu < mol.n_basis; mu++)
    {
        for (int nu = 0; nu < mol.n_basis; nu++)
        {
            z_a = mol.z_vals(mol.atom_map[mu]);
            z_b = mol.z_vals(mol.atom_map[nu]);
            x(mu, nu) = (mol.beta[z_a] + mol.beta[z_b]) * mol.density_matrix_total(mu, nu);
        }
    }

    return x;
}

//function to build the y matrix in the CNDO/2 analytical gradient expression
//inputs:
//      molecule mol - the molecule object with which to calculate the gradient
//outputs:
//      arma::mat - a matrix of dimension (n_atoms x n_atoms)
arma::mat gradient::build_y_matrix(molecule mol)
{
    arma::mat y(mol.n_atoms, mol.n_atoms);
    double sum;
  
    for (int atom_a = 0; atom_a < mol.n_atoms; atom_a++)
    {
        for (int atom_b = 0; atom_b < mol.n_atoms; atom_b++)
        {
            sum = 0;
            for (int mu = 0; mu < mol.n_basis; mu++)
            {
                for (int nu = 0; nu < mol.n_basis; nu++)
                {
                    if ((mol.atom_map[mu] == atom_a) && (mol.atom_map[nu] == atom_b))
                    {
                        sum += mol.density_matrix_alpha(mu, nu) * mol.density_matrix_alpha(mu, nu)
                                + mol.density_matrix_beta(mu, nu) * mol.density_matrix_beta(mu, nu);
                    }
                }
            }

            y(atom_a, atom_b) = mol.atomic_density(atom_a) * mol.atomic_density(atom_b)
                                - mol.valence_charge(atom_b) * mol.atomic_density(atom_a)
                                - mol.valence_charge(atom_a) * mol.atomic_density(atom_b)
                                - sum;

        }
    }

    return y;
}

//function to calculate the derivative of the nuclear repulsion energy wrt atomic perturbation
//inputs:
//      molecule mol - the molecule object with which to calculate the gradient
//outputs:
//      arma::mat - a matrix of dimension (3 x n_atoms)
arma::mat gradient::calculate_dv_dr(molecule mol)
{
    arma::mat dv_dr(mol.n_atoms, 3);
    arma::vec sum;
    arma::vec r;
    double r_norm;

    for (int atom_a = 0; atom_a < mol.n_atoms; atom_a++)
    {
        sum = {0, 0, 0};

        for (int atom_b = 0; atom_b < mol.n_atoms; atom_b++)
        {
            
            if (atom_a != atom_b)
            {
                r = mol.coords.row(atom_a).t() - mol.coords.row(atom_b).t();
                r_norm = arma::norm(r, 2);
                sum += - mol.valence_charge(atom_a) * mol.valence_charge(atom_b) * (1.0/2.0) * std::pow(r_norm, -3.0) * 2 *r * 27.211;
            }
            
            
        }

        dv_dr.row(atom_a) = sum.t();
    }

    return dv_dr;
}

//function to calculate the CNDO/2 analytical gradient expression
//inputs:
//      molecule mol - the molecule object with which to calculate the gradient
//outputs:
//      arma::mat - a matrix of dimension (3 x n_atoms)
arma::mat gradient::calculate_gradient(molecule mol)
{
    arma::mat x = gradient::build_x_matrix(mol);
    arma::mat y = gradient::build_y_matrix(mol);

    arma::field<arma::vec> dgamma = util::calculate_d_gamma_ab_matrix(mol);

    arma::field<arma::vec> ds = util::calculate_s_mu_nu_matrix(mol);

    arma::mat dv = gradient::calculate_dv_dr(mol);

    arma::vec sum = {0, 0, 0};
    arma::mat result(mol.n_atoms, 3);

    for (int atom_a = 0; atom_a < mol.n_atoms; atom_a++)
    {
        sum = {0, 0, 0};

        for (int mu = 0; mu < mol.n_basis; mu++)
        {
            for (int nu = 0; nu < mol.n_basis; nu++)
            {
                if ((mol.atom_map[mu] == atom_a) && (mol.atom_map[nu]!= atom_a))
                {
                    sum += x(mu, nu) * ds(mu, nu);
                }
            }
        }

        for (int atom_b = 0; atom_b < mol.n_atoms; atom_b++)
        {
            if (atom_a != atom_b)
            {
                sum += y(atom_a, atom_b) * dgamma(atom_a, atom_b);
            }
        }


        sum += dv.row(atom_a).t();

        result.row(atom_a) = sum.t();
    }

    return result;

}

//function to print the arma::field<arma::vec> objects to std::cout
//inputs:
//      arma::field<arma::vec> field - the given field
//outputs: none. print to std::cout
void gradient::print_field(arma::field<arma::vec> field)
{

    for (int atom_a = 0; atom_a < field.n_rows; atom_a++)
    {
        for (int atom_b = 0; atom_b < field.n_cols; atom_b++)
        {
            std::cout << field(atom_a, atom_b) << "------" << std::endl;
        }
    }
}