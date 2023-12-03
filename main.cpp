#include"integrator.h"
#include"gradient.h"
#include"util.h"
#include"SCF.h"
#include"atomic_orbital.h"
#include"molecule.h"
#include"primitive_guassian.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>


int main() {

    molecule mol_c2h4("Molecules/C2H4.txt", 6, 6);
    std::cout << "Molecule: C2H4" << std::endl;
    //mol_h2.print_info();
    SCF::DIIS(mol_c2h4, 10e-6, 12, true);
    
    arma::mat grad = gradient::calculate_gradient(mol_c2h4);
    std::cout << "Gradient" << std::endl;
    std::cout << grad << std::endl;


    return 0;
}
