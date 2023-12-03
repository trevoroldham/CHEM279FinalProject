#include"util.h"
#include"SCF.h"
#include"gradient.h"
#include"molecule.h"
#include"integrator.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>


//function to perform the fixed point iterations of the self consistent field for CNDO/2
//inputs:
//      bool verbose - if TRUE, prints information about current iteration. If false, no print output
//outputs: none
void SCF::fixed_point_iteration(molecule & mol, double tolerance, bool verbose)
{
    int iter = 0;

    arma::mat p_alpha_old;
    p_alpha_old.resize(mol.n_basis, mol.n_basis);
    arma::mat p_beta_old;
    p_beta_old.resize(mol.n_basis, mol.n_basis);

    mol.coefficient_matrix_alpha = mol.calculate_coefficient_matrix(mol.fock_matrix_alpha);
    mol.coefficient_matrix_beta = mol.calculate_coefficient_matrix(mol.fock_matrix_beta);

    if (verbose)
    {
        std::cout << "Performing SCF Iterations (Fixed Point)" << std::endl;
    }

    do {
        if (verbose)
        {
            std::cout << "Iteration: " << iter << std::endl;
        }
        arma::vec eigenvalues;
        arma::mat eigenvectors;



        p_alpha_old = mol.density_matrix_alpha;
        mol.density_matrix_alpha = mol.calculate_density_matrix(mol.coefficient_matrix_alpha, mol.p_electrons);

        p_beta_old = mol.density_matrix_beta;

        if (mol.q_electrons > 0)
        {
            mol.density_matrix_beta = mol.calculate_density_matrix(mol.coefficient_matrix_beta, mol.q_electrons);
        }

        mol.density_matrix_total = mol.calculate_total_density_matrix(mol.density_matrix_alpha, mol.density_matrix_beta, mol.atomic_density);


        mol.fock_matrix_alpha = mol.calculate_fock_matrix(mol.density_matrix_alpha);
        arma::eig_sym(eigenvalues, eigenvectors, mol.fock_matrix_alpha);
        mol.coefficient_matrix_alpha = eigenvectors;

        
        if (mol.q_electrons > 0)
        {
            

            mol.fock_matrix_beta = mol.calculate_fock_matrix(mol.density_matrix_beta);
            arma::eig_sym(eigenvalues, eigenvectors, mol.fock_matrix_beta);
            mol.coefficient_matrix_beta = eigenvectors;
        }
        
        else 
        {
            mol.density_matrix_total = mol.density_matrix_alpha;

        }

        iter++;
       

    } while ((arma::norm(mol.density_matrix_alpha - p_alpha_old, 2) > tolerance) || (arma::norm(mol.density_matrix_beta - p_beta_old, 2) > tolerance));

}

void SCF::DIIS(molecule & mol, double tolerance, int history, bool verbose)
{
    int iter = 0;

    std::vector<arma::mat> fock_matrices_alpha(history);
    std::vector<arma::mat> fock_matrices_beta(history);
    std::vector<arma::mat> error_matrices_alpha(history);
    std::vector<arma::mat> error_matrices_beta(history);

    arma::mat p_alpha_old;
    p_alpha_old.resize(mol.n_basis, mol.n_basis);
    arma::mat p_beta_old;
    p_beta_old.resize(mol.n_basis, mol.n_basis);

    mol.coefficient_matrix_alpha = mol.calculate_coefficient_matrix(mol.fock_matrix_alpha);
    mol.coefficient_matrix_beta = mol.calculate_coefficient_matrix(mol.fock_matrix_beta);

    if (verbose)
    {
        std::cout << "Performing SCF Iterations (DIIS)" << std::endl;
    }

    do {
        if (verbose)
        {
            std::cout << "Iteration: " << iter << std::endl;
        }
        arma::vec eigenvalues;
        arma::mat eigenvectors;



        p_alpha_old = mol.density_matrix_alpha;
        mol.density_matrix_alpha = mol.calculate_density_matrix(mol.coefficient_matrix_alpha, mol.p_electrons);

        p_beta_old = mol.density_matrix_beta;

        if (mol.q_electrons > 0)
        {
            mol.density_matrix_beta = mol.calculate_density_matrix(mol.coefficient_matrix_beta, mol.q_electrons);
        }

        mol.density_matrix_total = mol.calculate_total_density_matrix(mol.density_matrix_alpha, mol.density_matrix_beta, mol.atomic_density);


        mol.fock_matrix_alpha = mol.calculate_fock_matrix(mol.density_matrix_alpha);
        fock_matrices_alpha[iter] = mol.fock_matrix_alpha;
        error_matrices_alpha[iter] = mol.fock_matrix_alpha * mol.density_matrix_alpha - mol.density_matrix_alpha*mol.fock_matrix_alpha;

        arma::eig_sym(eigenvalues, eigenvectors, mol.fock_matrix_alpha);
        mol.coefficient_matrix_alpha = eigenvectors;

        
        if (mol.q_electrons > 0)
        {
            

            mol.fock_matrix_beta = mol.calculate_fock_matrix(mol.density_matrix_beta);
            fock_matrices_beta[iter] = mol.fock_matrix_beta;
            error_matrices_beta[iter] = mol.fock_matrix_beta * mol.density_matrix_beta - mol.density_matrix_beta*mol.fock_matrix_beta;

            arma::eig_sym(eigenvalues, eigenvectors, mol.fock_matrix_beta);
            mol.coefficient_matrix_beta = eigenvectors;
        }
        
        else 
        {
            mol.density_matrix_total = mol.density_matrix_alpha;
            fock_matrices_beta[iter] = mol.fock_matrix_beta;
            error_matrices_beta[iter] = mol.fock_matrix_beta * mol.density_matrix_beta - mol.density_matrix_beta*mol.fock_matrix_beta;

        }

        iter++;
       

    } while (((arma::norm(mol.density_matrix_alpha - p_alpha_old, 2) > tolerance) || (arma::norm(mol.density_matrix_beta - p_beta_old, 2) > tolerance)) && (iter < history));

    arma::mat grad = gradient::calculate_gradient(mol);
    std::cout << "Gradient" << std::endl;
    std::cout << grad << std::endl;
    //break if tolerance has already been reached
    if (iter < history)
    {
        return;
    }

    //begin DIIS iterations
    if (verbose)
    {
        std::cout << "Beginning DIIS iterations" << std::endl;
    }

    //initialize field to hold fock matrix  and error matrix flattened
    std::vector<arma::vec> f_alpha(history);
    std::vector<arma::vec> e_alpha(history);
    std::vector<arma::vec> f_beta(history);
    std::vector<arma::vec> e_beta(history);

    //initialize DIIS matrices
    arma::mat DIIS_matrix_alpha(history + 1, history + 1);
    arma::mat DIIS_matrix_beta(history + 1, history + 1);

    do {

        if (verbose)
        {
            std::cout << "Iteration: " << iter << std::endl;
        }

        //save old density matrices
        p_alpha_old = mol.density_matrix_alpha;
        p_beta_old = mol.density_matrix_beta;
        mol.density_matrix_alpha = mol.calculate_density_matrix(mol.coefficient_matrix_alpha, mol.p_electrons);

        //set vector elements of f_alpha and e_alpha
        for (int i = 0; i < history; i++)
        {
            f_alpha[i] = arma::vectorise(fock_matrices_alpha[i]);
            e_alpha[i] = arma::vectorise(error_matrices_alpha[i]);
        }

        //set values of DIIS matrix
        for (int i = 0; i < history + 1; i++)
        {
            for (int j = 0; j < history + 1; j++)
            {
                if (i == history || j == history)
                {
                    DIIS_matrix_alpha(i, j) = -1;
                }

                if (i == history && j == history)
                {
                    DIIS_matrix_alpha(i, j) = 0;
                }

                if (i != history && j != history)
                {
                    DIIS_matrix_alpha(i, j) = arma::dot(e_alpha[i], e_alpha[j]);
                }
           
            }
        }

        
        //initialize the y vector in the DIIS matrix equation and solve for coefficients
        arma::vec y(history + 1);
        y(history) = -1;
        arma::vec sol = arma::pinv(DIIS_matrix_alpha) * y;
        arma::mat f_star(mol.n_basis, mol.n_basis);

        //compute f_star using coefficients and fock_matrices_alpha
        for (int i = 0; i < history; i++)
        {
            f_star += sol(i) * fock_matrices_alpha[i];
        }

        //solve the symmetric eigenvalue problem for the molecular coefficients
        arma::vec eigenvalues;
        arma::mat eigenvectors;
        arma::eig_sym(eigenvalues, eigenvectors, f_star);
        mol.coefficient_matrix_alpha = eigenvectors;

        //save density matrix and fock matrix in molecule
        mol.fock_matrix_alpha = f_star;
        
        //shift the entries in past fock history down one index
        for (int i = 0; i < history - 1; i++)
        {
            fock_matrices_alpha[i] = fock_matrices_alpha[i+1];
            error_matrices_alpha[i] = error_matrices_alpha[i+1];
        }

        //set the new most recent fock matrix as f_star
        fock_matrices_alpha[history - 1] = f_star;
        error_matrices_alpha[history - 1] = mol.fock_matrix_alpha * mol.density_matrix_alpha - mol.density_matrix_alpha * mol.fock_matrix_alpha;


        if (mol.q_electrons > 0)
        {
            mol.density_matrix_beta = mol.calculate_density_matrix(mol.coefficient_matrix_beta, mol.q_electrons);

            //set vector elements of f_alpha and e_alpha
            for (int i = 0; i < history; i++)
            {
                f_beta[i] = arma::vectorise(fock_matrices_beta[i]);
                e_beta[i] = arma::vectorise(error_matrices_beta[i]);
            }

            //set values of DIIS matrix
            for (int i = 0; i < history + 1; i++)
            {
                for (int j = 0; j < history + 1; j++)
                {
                    if (i == history || j == history)
                    {
                        DIIS_matrix_beta(i, j) = -1;
                    }

                    if (i == history && j == history)
                    {
                        DIIS_matrix_beta(i, j) = 0;
                    }

                    if (i != history && j != history)
                    {
                        DIIS_matrix_beta(i, j) = arma::dot(e_beta[i], e_beta[j]);
                    }
            
                }
            }

            //initialize the y vector in the DIIS matrix equation and solve for coefficients
            arma::vec y(history + 1);
            y(history) = -1;
            sol = arma::pinv(DIIS_matrix_beta) * y;
            arma::mat f_star(mol.n_basis, mol.n_basis);

            //compute f_star using coefficients and fock_matrices_beta
            for (int i = 0; i < history; i++)
            {
                f_star += sol(i) * fock_matrices_beta[i];
            }

            //solve the symmetric eigenvalue problem for the molecular coefficients
            arma::vec eigenvalues;
            arma::mat eigenvectors;
            arma::eig_sym(eigenvalues, eigenvectors, f_star);
            mol.coefficient_matrix_beta = eigenvectors;

            //save density matrix and fock matrix in molecule
            mol.fock_matrix_beta = f_star;
            
            //shift the entries in past fock history down one index
            for (int i = 0; i < history - 1; i++)
            {
                fock_matrices_beta[i] = fock_matrices_beta[i+1];
                error_matrices_beta[i] = error_matrices_beta[i+1];
            }

            //set the new most recent fock matrix as f_star
            fock_matrices_beta[history - 1] = f_star;
            error_matrices_beta[history - 1] = mol.fock_matrix_beta * mol.density_matrix_beta - mol.density_matrix_beta * mol.fock_matrix_beta;

        }

        if (mol.q_electrons < 1) 
        {
            mol.density_matrix_total = mol.density_matrix_alpha;

        }
        else
        {
            mol.density_matrix_total = mol.calculate_total_density_matrix(mol.density_matrix_alpha, mol.density_matrix_beta, mol.atomic_density);

        }

        iter++;
        //std::cout << "norm" << std::endl;
        //std::cout << arma::norm(mol.density_matrix_alpha - p_alpha_old, 2) << std::endl;
        //std::cout << arma::norm(mol.density_matrix_beta - p_beta_old, 2) << std::endl;

    } while ((arma::norm(mol.density_matrix_alpha - p_alpha_old, 2) > tolerance) || (arma::norm(mol.density_matrix_beta - p_beta_old, 2) > tolerance));
 
}
