#include"util.h"
#include"SCF.h"
#include"gradient.h"
#include"molecule.h"
#include"integrator.h"
#include"BFGS.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>

//function to perform BFGS minimization of the gradient
//inputs:
//      molecule mol - the input molecule
//      double tolerance - the desired tolerance for the end condition
//      int max_iter - the maximum number of iterations
void BFGS::minimize(molecule & mol, double tolerance, int max_iter)
{
    arma::vec x_0 = arma::vectorise(mol.coords);
    int size = x_0.size();
    int iter = 0;

    //set vectors s and y 
    arma::mat s(mol.n_atoms * 3, 1);
    arma::mat y(mol.n_atoms * 3, 1);

    //initialize hessian to the identity matrix
    std::cout << "initializing hessian" << std::endl;
    arma::mat hess(size, size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                hess(i, j) = 1;
            }

            else
            {
                hess(i, j) = 0;
            }
        }
    }

    std::cout << "initializing p" << std::endl;
    arma::vec p = -arma::vectorise(gradient::calculate_gradient(mol));

    double alpha = 1;
    arma::mat x = x_0;
    std::cout << "x_0" << std::endl;
    std::cout << x << std::endl;

    int k = 0;
    std::cout << "initializing grad_old" << std::endl;
    arma::vec grad_old = arma::vectorise(gradient::calculate_gradient(mol));
    std::cout << "Grad old" << std::endl;
    std::cout << grad_old << std::endl;
    
    arma::vec grad_new(mol.n_atoms * 3);

    while ((arma::norm(gradient::calculate_gradient(mol), 2) > tolerance) & (iter < max_iter))
    {
        //x = x.resize(mol.n_atoms * 3, 1);
        std::cout << "making grad_old" << std::endl;
        grad_old = arma::vectorise(gradient::calculate_gradient(mol));
        std::cout << "making p" << std::endl;
        p = arma::pinv(hess) * -grad_old;

        //find alpha with line search
        std::cout << "Initializing a" << std::endl;
        alpha = BFGS::line_search(mol, p, grad_old, 100);
        s = alpha*p;
        std::cout << "S" << std::endl;
        std::cout << s << std::endl;
        std::cout << "x +  s" << std::endl;
        std::cout << x + s << std::endl;

        arma::mat x_new = x + s;
        x_new = arma::reshape(x_new, mol.n_atoms, 3);
        std::cout << "x_new" << std::endl;
        std::cout << x_new << std::endl;
        mol.coords = x_new;
        SCF::fixed_point_iteration(mol, 1e-5, false);

        grad_new = arma::vectorise(gradient::calculate_gradient(mol));

        y = grad_new - grad_old;
        std::cout << "grad_new" << std::endl;
        std::cout << grad_new << std::endl;
        std::cout << "y" << std::endl;
        std::cout << y << std::endl;
        std::cout << "summing hessian" << std::endl;
        //hess = hess + (y * y.t()) / arma::dot(y.t(), s) - (hess *  s * s.t() * hess.t()) / (s.t() * hess * s);
        //hess = hess + (arma::dot(y, y.t()) / (arma::dot(y.t(), s))) - (hess * arma::dot(s, s.t()) * hess.t()) / arma::as_scalar(s.t() * hess * s);
        arma::mat ys = y * s.t();
        arma::mat Hs = hess * s;
        //arma::mat sHy = s.t() * Hy;
        hess = hess + (y * y.t()) / arma::dot(y.t(), s) - (Hs * Hs.t()) / arma::dot(s.t(), Hs);

        std::cout << "Hessian " << iter << std::endl;
        std::cout << hess << std::endl;
        x = x_new;
        x = arma::vectorise(x);
        iter++;

    }

    //x.resize(mol.n_atoms, 3);
    return;

    //repeat
}


//function to compute optimal value of alpha for a line search along direction of the gradient
//inputs:
//      molecule mol - the input molecule
//      arma::vec p - the direction of the gradient
//      arma::vec grad - the gradient
//outputs:
//      double a - the value of alpha step along the direction of the gradient
double BFGS::line_search(molecule mol, arma::vec p, arma::vec grad, int max_iter)
{
    int iter = 0;
    double a = 10;
    double c1 = 1e-4;
    double c2 = 0.9;

    std::cout << "Initializing x_0" << std::endl;
    arma::vec x_0 = arma::vectorise(mol.coords);
    double fx = mol.calculate_energy(false);

    std::cout << "Initializing x_new" << std::endl;
    arma::mat x_new = x_0 + a * p;
    std::cout << "resizing x_new" << std::endl;
    x_new = arma::reshape(x_new, mol.n_atoms, 3);
    std::cout << "setting mol.coords" << std::endl;
    mol.coords = x_new;
    SCF::fixed_point_iteration(mol, 1e-5, false);

    std::cout << "Initializing grad_new" << std::endl;
    arma::vec grad_new = arma::vectorise(gradient::calculate_gradient(mol));


    while (((mol.calculate_energy(false)) >= (fx + c1 * a * arma::dot(grad.t(), p))) || (arma::dot(grad_new.t(), p) <= c2 * arma::dot(grad.t(), p)))
    {
        a = a * 0.9;
        x_new = arma::vectorise(x_new);
        x_new = x_0 + a * p;
        x_new = arma::reshape(x_new, mol.n_atoms, 3);
        mol.coords = x_new;
        SCF::fixed_point_iteration(mol, 1e-5, false);
        grad_new = arma::vectorise(gradient::calculate_gradient(mol));

        iter++;
        if (iter >= max_iter)
        {
            break;
        }
    }

    std::cout << "Iter" << iter << std::endl;
    std::cout << "alpha = " << a << std::endl;
    return a;

}

