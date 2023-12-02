#ifndef MOLECULE_H
#define MOLECULE_H
#include"primitive_guassian.h"
#include"atomic_orbital.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>
#include<fstream>

class molecule {

    public:
        //the number of atoms in molecule
        int n_atoms;

        //the x,y,z coordinates of the atoms in molecule
        arma::mat coords;

        //the integer values of nuclear atomic number for each atom in molecule
        arma::ivec z_vals;

        //the integer value of the valence atomic number for each atom
        arma::ivec valence_charge;

        //the number of basis functions in the molecule
        int n_basis;

        //number of p and q electrons (passed as input to molecule constructor)
        int n_electrons;
        int p_electrons;
        int q_electrons;

        //the molecular charge (passed as input from file of coordinates)
        int charge;

        //the overlap matrix calculated using functions in util and integrator
        arma::mat overlap_matrix;

        //the matrices used in the CNDO/2 calculation
        arma::mat fock_matrix_alpha;
        arma::mat fock_matrix_beta;
        arma::mat coefficient_matrix_alpha;
        arma::mat coefficient_matrix_beta;
        arma::mat density_matrix_alpha;
        arma::mat density_matrix_beta;
        arma::mat density_matrix_total;
        arma::mat gamma_matrix;
        arma::mat core_hamiltonian;

        //a vector containing the total electron density on each atom
        arma::vec atomic_density;

        //maps containing the CNDO/2 parameters for each orbital type
        std::map<int, double> s_parameters;
        std::map<int, double> p_parameters;

        //map containing the beta values for each atom type
        std::map<int, double> beta;

        //map from atomic char symbols to integer atomic numbers
        std::map<char, int> z_symbols;

        //functions to calculate the matrices in the CNDO/2 model
        arma::mat calculate_fock_matrix(arma::mat density_matrix);
        arma::mat calculate_core_hamiltonian();
        arma::mat calculate_coefficient_matrix(arma::mat fock_matrix);
        arma::mat calculate_density_matrix(arma::mat coeff_matrix, int occ_mos);
        arma::mat calculate_total_density_matrix(arma::mat density_alpha, arma::mat density_beta, arma::vec & atomic_density);

        //function to iterate over the self consistent field method through fixed point iteration
        void scf_steps(bool verbose);

        //function to calculate energy in electron volts after performing SCF minimization
        double calculate_energy(bool verbose);

        //function to print info about the matrix objects in molecule
        void print_info();

        //vector containing the atomic orbitals in the molecule
        std::vector<atomic_orbital> orbitals;

        //vector containing the integer atomic number of a given orbital
        std::vector<int> orbital_species;

        //map containing the associated atom index of a given atomic orbital
        std::map<int, int> atom_map;

        //map containing the char value (s or p) to a given integer index of an atomic orbital
        std::map<int, char> orbital_type;

        //default constructor
        molecule(const char *filename, int p, int q);
};

#endif