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
        int n_atoms;
        arma::mat coords;
        arma::ivec z_vals;
        arma::ivec valence_charge;
        int n_basis;
        int n_electrons;
        int p_electrons;
        int q_electrons;
        int charge;
        arma::mat overlap_matrix;
        arma::mat fock_matrix_alpha;
        arma::mat fock_matrix_beta;
        arma::mat coefficient_matrix_alpha;
        arma::mat coefficient_matrix_beta;
        arma::mat density_matrix_alpha;
        arma::mat density_matrix_beta;
        arma::mat density_matrix_total;
        arma::vec atomic_density;
        arma::mat gamma_matrix;
        arma::mat core_hamiltonian;

        std::map<int, double> s_parameters;
        std::map<int, double> p_parameters;
        std::map<int, double> beta;
        std::map<char, int> z_symbols;

        arma::mat calculate_fock_matrix(arma::mat density_matrix);
        arma::mat calculate_core_hamiltonian();
        arma::mat calculate_coefficient_matrix(arma::mat fock_matrix);
        arma::mat calculate_density_matrix(arma::mat coeff_matrix, int occ_mos);
        arma::mat calculate_total_density_matrix(arma::mat density_alpha, arma::mat density_beta, arma::vec & atomic_density);
        void scf_steps(bool verbose);
        double calculate_energy(bool verbose);
        void print_info();

        std::vector<atomic_orbital> orbitals;
        std::vector<int> orbital_species;
        std::map<int, int> atom_map;
        std::map<int, char> orbital_type;

        molecule(const char *filename, int p, int q);
};

#endif