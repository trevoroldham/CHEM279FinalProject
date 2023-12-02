#include"integrator.h"
#include"gradient.h"
#include"util.h"
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
    mol_c2h4.scf_steps(false);
    
    arma::mat grad = gradient::calculate_gradient(mol_c2h4);
    std::cout << "Gradient" << std::endl;
    std::cout << grad << std::endl;

    molecule mol_hf("Molecules/HF.txt", 4, 3);
    std::cout << "Molecule: HF" << std::endl;
    //mol_h2.print_info();
    mol_hf.scf_steps(false);
    
    grad = gradient::calculate_gradient(mol_hf);
    std::cout << "Gradient" << std::endl;
    std::cout << grad << std::endl;

    return 0;
}
