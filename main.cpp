#include"integrator.h"
#include"gradient.h"
#include"util.h"
#include"SCF.h"
#include"atomic_orbital.h"
#include"molecule.h"
#include"BFGS.h"
#include"primitive_guassian.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>


int main() {

    molecule mol_c2h4("Molecules/C2H4.txt", 6, 6);
    std::cout << "Molecule: C2H4" << std::endl;
    //mol_h2.print_info();
    SCF::fixed_point_iteration(mol_c2h4, 1e-6, true);
    
    std::cout << "Energy: " << mol_c2h4.calculate_energy(false) << std::endl;
    arma::mat grad = gradient::calculate_gradient(mol_c2h4);
    std::cout << "Gradient" << std::endl;
    std::cout << grad << std::endl;

    BFGS::minimize(mol_c2h4, 1e-5, 10);
    std::cout << "Coordinates after minimization" << std::endl;
    std::cout << mol_c2h4.coords << std::endl;
    std::cout << "Gradient after minimization" << std::endl;
    SCF::fixed_point_iteration(mol_c2h4, 1e-6, false);
    grad = gradient::calculate_gradient(mol_c2h4);
    std::cout << grad << std::endl;
    std::cout << "Energy after minimization" << std::endl;
    std::cout << mol_c2h4.calculate_energy(false) << std::endl;


    return 0;
}
