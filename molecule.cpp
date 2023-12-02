#include"primitive_guassian.h"
#include"molecule.h"
#include"util.h"
#include<iostream>
#include<armadillo>
#include<cstdio>
#include<cmath>
#include<vector>
#include<fstream>
#include<cassert>

//constructor reads from txt file and initializes relevant member variables
molecule::molecule(const char *filename, int p, int q)
{
    p_electrons = p;
    q_electrons = q;
    z_symbols['H'] = 1;
    z_symbols['C'] = 6;
    z_symbols['N'] = 7;
    z_symbols['O'] = 8;
    z_symbols['F'] = 9;

    
    std::ifstream input_stream(filename);

    assert(input_stream.good());

    input_stream >> n_atoms;

    char symbols[n_atoms];

    input_stream >> charge;

    z_vals.resize(n_atoms);
    coords.resize(n_atoms, 3);

    for (int i = 0; i < n_atoms; i++)
    {
        input_stream >> symbols[i] >> coords(i, 0) >> coords(i, 1) >> coords(i, 2);
    }

    input_stream.close();

    for (int i = 0; i < n_atoms; i++)
    {
        z_vals(i) = z_symbols[symbols[i]];
    }

    valence_charge.resize(n_atoms);
    for (int i = 0; i < n_atoms; i++)
    {
        if (z_vals(i) != 1)
        {
            valence_charge(i) = z_vals(i) - 2;
        }

        else
        {
            valence_charge(i) = z_vals(i);
        }
    }

    //define coefficients and exponents for orbital guassians
    arma::vec hydrogen_coefficients(3);
    hydrogen_coefficients = {0.15432897, 0.53532814, 0.44463454};
    arma::vec hydrogen_exponents(3);
    hydrogen_exponents = {3.42525091, 0.62391373, 0.16885540};

    arma::vec carbon_exponents(3);
    carbon_exponents = {2.94124940, 0.68348310, 0.22228990};
    arma::vec carbon_2s_coefficients(3);
    carbon_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
    arma::vec carbon_2p_coefficients(3);
    carbon_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};

    arma::vec nitrogen_exponents(3);
    nitrogen_exponents = {3.78045590, 0.87849660, 0.28571440};
    arma::vec nitrogen_2s_coefficients(3);
    nitrogen_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
    arma::vec nitrogen_2p_coefficients(3);
    nitrogen_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};

    arma::vec oxygen_exponents(3);
    oxygen_exponents = {5.03315130, 1.16959610, 0.38038900};
    arma::vec oxygen_2s_coefficients(3);
    oxygen_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
    arma::vec oxygen_2p_coefficients(3);
    oxygen_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};

    arma::vec flourine_exponents(3);
    flourine_exponents = {6.46480320, 1.50228120, 0.48858850};
    arma::vec flourine_2s_coefficients(3);
    flourine_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
    arma::vec flourine_2p_coefficients(3);
    flourine_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};



    //count the number of each type of atom
    int n_electrons;
    int carbon_count = 0;
    int hydrogen_count = 0;
    int nitrogen_count = 0;
    int oxygen_count = 0;
    int flourine_count = 0;

    //count the current number of orbitals
    int orbital_count = 0;

    //add the orbitals of each atom to the molecule
    //add the type of atom to orbital_species vector
    //add the type of atom for each orbital to atom_map
    for (int i = 0; i < n_atoms; i++)
    {
        if (z_vals[i] == 6)
        {
            carbon_count++;
            n_electrons += 4;
            atomic_orbital orbital_1 = atomic_orbital(coords.row(i).t(), 0, 0, 0, carbon_exponents, carbon_2s_coefficients);
            atomic_orbital orbital_2 = atomic_orbital(coords.row(i).t(), 1, 0, 0, carbon_exponents, carbon_2p_coefficients);
            atomic_orbital orbital_3 = atomic_orbital(coords.row(i).t(), 0, 1, 0, carbon_exponents, carbon_2p_coefficients);
            atomic_orbital orbital_4 = atomic_orbital(coords.row(i).t(), 0, 0, 1, carbon_exponents, carbon_2p_coefficients);

            orbitals.push_back(orbital_1);
            orbital_species.push_back(6);
            orbital_type[orbital_count] = 'S';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_2);
            orbital_species.push_back(6);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_3);
            orbital_species.push_back(6);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_4);
            orbital_species.push_back(6);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
        }

        else if (z_vals[i] == 1)
        {
            hydrogen_count++;
            n_electrons += 1;
            atomic_orbital orbital_1 = atomic_orbital(coords.row(i).t(), 0, 0, 0, hydrogen_exponents, hydrogen_coefficients);
            orbitals.push_back(orbital_1);
            orbital_species.push_back(1);
            orbital_type[orbital_count] = 'S';
            atom_map[orbital_count] = i;
            orbital_count++;
        }

        else if (z_vals[i] == 7)
        {
            nitrogen_count++;
            n_electrons += 5;
            atomic_orbital orbital_1 = atomic_orbital(coords.row(i).t(), 0, 0, 0, nitrogen_exponents, nitrogen_2s_coefficients);
            atomic_orbital orbital_2 = atomic_orbital(coords.row(i).t(), 1, 0, 0, nitrogen_exponents, nitrogen_2p_coefficients);
            atomic_orbital orbital_3 = atomic_orbital(coords.row(i).t(), 0, 1, 0, nitrogen_exponents, nitrogen_2p_coefficients);
            atomic_orbital orbital_4 = atomic_orbital(coords.row(i).t(), 0, 0, 1, nitrogen_exponents, nitrogen_2p_coefficients);

            orbitals.push_back(orbital_1);
            orbital_species.push_back(7);
            orbital_type[orbital_count] = 'S';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_2);
            orbital_species.push_back(7);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_3);
            orbital_species.push_back(7);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_4);
            orbital_species.push_back(7);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
        }

        else if (z_vals[i] == 8)
        {
            oxygen_count++;
            n_electrons += 6;
            atomic_orbital orbital_1 = atomic_orbital(coords.row(i).t(), 0, 0, 0, oxygen_exponents, oxygen_2s_coefficients);
            atomic_orbital orbital_2 = atomic_orbital(coords.row(i).t(), 1, 0, 0, oxygen_exponents, oxygen_2p_coefficients);
            atomic_orbital orbital_3 = atomic_orbital(coords.row(i).t(), 0, 1, 0, oxygen_exponents, oxygen_2p_coefficients);
            atomic_orbital orbital_4 = atomic_orbital(coords.row(i).t(), 0, 0, 1, oxygen_exponents, oxygen_2p_coefficients);

            orbitals.push_back(orbital_1);
            orbital_species.push_back(8);
            orbital_type[orbital_count] = 'S';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_2);
            orbital_species.push_back(8);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_3);
            orbital_species.push_back(8);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_4);
            orbital_species.push_back(8);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
        }

        else if (z_vals[i] == 9)
        {
            flourine_count++;
            n_electrons += 7;
            atomic_orbital orbital_1 = atomic_orbital(coords.row(i).t(), 0, 0, 0, flourine_exponents, flourine_2s_coefficients);
            atomic_orbital orbital_2 = atomic_orbital(coords.row(i).t(), 1, 0, 0, flourine_exponents, flourine_2p_coefficients);
            atomic_orbital orbital_3 = atomic_orbital(coords.row(i).t(), 0, 1, 0, flourine_exponents, flourine_2p_coefficients);
            atomic_orbital orbital_4 = atomic_orbital(coords.row(i).t(), 0, 0, 1, flourine_exponents, flourine_2p_coefficients);

            orbitals.push_back(orbital_1);
            orbital_species.push_back(9);
            orbital_type[orbital_count] = 'S';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_2);
            orbital_species.push_back(9);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_3);
            orbital_species.push_back(9);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
            orbitals.push_back(orbital_4);
            orbital_species.push_back(9);
            orbital_type[orbital_count] = 'P';
            atom_map[orbital_count] = i;
            orbital_count++;
        }
    }

    //adjust n_electrons using value of molecular charge
    n_electrons = n_electrons - charge;


    //compute the number of basis functions in the molecule
    n_basis = 4*carbon_count + 4*oxygen_count + 4*nitrogen_count + 4*flourine_count +  hydrogen_count;

    //compute normalization constants for orbital guassians

    
    for (int i = 0; i < n_basis; i++)
    {
        orbitals[i].set_normalization_constants();
    }


    //compute the overlap matrix
    overlap_matrix.resize(n_basis, n_basis);

    for (int i = 0; i < n_basis; i++)
    {
        for (int j = 0; j < n_basis; j++)
        {
            overlap_matrix(i, j) = util::calculate_overlap(orbitals[i], orbitals[j]);
        }
    }

    

    //compute the gamma matrix (A x A) using only s orbitals from each atom
    gamma_matrix.resize(n_atoms, n_atoms);
    int atom_a = 0;
    int orbital_a_index = 0;
    int atom_b = 0;
    int orbital_b_index = 0;
    double gamma = 0;

    while (atom_a < n_atoms)
    {
        atom_b  = 0;
        orbital_b_index = 0;

        while (atom_b < n_atoms)
        {
            atomic_orbital a = orbitals[orbital_a_index];
            atomic_orbital b = orbitals[orbital_b_index];

            gamma = util::calculate_gamma(a, b);
            gamma_matrix(atom_a, atom_b) = gamma;

            if (z_vals(atom_b) != 1)
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

        if (z_vals(atom_a) != 1)
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

    


    //create the density matrices and initialize them to the zero matrix
    density_matrix_alpha.resize(n_basis, n_basis);
    density_matrix_beta.resize(n_basis, n_basis);
    density_matrix_total.resize(n_basis, n_basis);

    for (int i = 0; i < n_basis; i++)
    {
        for (int j = 0; j < n_basis; j++)
        {
            density_matrix_alpha(i, j) = 0;
        }
        
    }

    for (int i = 0; i < n_basis; i++)
    {
        for (int j = 0; j < n_basis; j++)
        {
            density_matrix_beta(i, j) = 0;
        }
        
    }

    for (int i = 0; i < n_basis; i++)
    {
        for (int j = 0; j < n_basis; j++)
        {
            density_matrix_total(i, j) = 0;
        }
        
    }

    
    //calculate density for each atom and save as a vector
    atomic_density.resize(n_atoms);
    for (int mu = 0; mu < n_basis; mu++)
    {
        atomic_density(atom_map[mu]) += density_matrix_total(mu, mu);
    }

    
    //define paramters (1/2)(I_s + A_s) and (1/2)(I_p + A_p)
    s_parameters[1] = 7.176;
    s_parameters[6] = 14.051;
    s_parameters[7] = 19.316;
    s_parameters[8] = 25.390;
    s_parameters[9] = 32.272;

    p_parameters[6] = 5.572;
    p_parameters[7] = 7.275;
    p_parameters[8] = 9.111;
    p_parameters[9] = 11.080;

    beta[1] = -9;
    beta[6] = -21;
    beta[7] = -25;
    beta[8] = -31;
    beta[9] = -39;

    //calculate initial matrices before the SCF steps
    fock_matrix_alpha = calculate_fock_matrix(density_matrix_alpha);
    fock_matrix_beta = calculate_fock_matrix(density_matrix_beta);
    
    core_hamiltonian = calculate_core_hamiltonian();
}

//define functions to calculate fock matrix, density matrix, and coefficient matrix
arma::mat molecule::calculate_fock_matrix(arma::mat density_matrix)
{
    arma::mat matrix;
    double constant;
    double b;
    char type_mu;
    char type_nu;
    double sum;
    int z_a;
    int z_b;

    matrix.resize(n_basis, n_basis);

    for (int mu = 0; mu < n_basis; mu++)
    {
        //get valence charge of nucleus
        if (z_vals(atom_map[mu]) != 1)
        {
            z_a = z_vals(atom_map[mu]) - 2;
        }

        else
        {
            z_a = 1;    
        }

        //get type of orbital mu (s or p)
        type_mu = orbital_type[mu];

        for (int nu = 0; nu < n_basis; nu++)
        {
            //diagonal elements
            
            if (mu == nu)
            {
                //get paramters for orbital mu, s or p
                if (type_mu == 'S')
                {
                    constant = - s_parameters[orbital_species[mu]];
                }
                else if (type_mu == 'P')
                {
                    constant = - p_parameters[orbital_species[mu]];
                }

                //calculate sum in the expression for diagonal fock matrix elements
                sum=0;

               for (int atom_b = 0; atom_b < n_atoms; atom_b++)
               {
                    //get valence charge of nucleus
                    if (z_vals(atom_b) != 1)
                    {
                        z_b = z_vals(atom_b) - 2;
                    }

                    else 
                    {
                        z_b = 1;
                    
                    }

                    if (atom_map[mu] != atom_b)
                    {
                        sum += (atomic_density(atom_b) - z_b) * gamma_matrix(atom_map[mu], atom_b);
                    }
               }

                matrix(mu, nu) = constant + 
                                            (
                                                (atomic_density(atom_map[mu]) - z_a) - 
                                                (density_matrix(mu, mu) - (1.0/2.0))
                                            ) * 
                                                gamma_matrix(atom_map[mu], atom_map[mu]) + sum;
                                            
            }

            else 
            {
                matrix(mu, nu) = (1.0/2.0) * (beta[orbital_species[mu]] + beta[orbital_species[nu]])
                                *overlap_matrix(mu, nu) - 
                                density_matrix(mu, nu) * gamma_matrix(atom_map[mu], atom_map[nu]);
                                
            }
        }
    }

    return matrix;
}

arma::mat molecule::calculate_core_hamiltonian()
{
    arma::mat matrix;
    double constant;
    double b;
    char type_mu;
    char type_nu;
    double sum;
    int z_a;
    int z_b;

    matrix.resize(n_basis, n_basis);

    for (int mu = 0; mu < n_basis; mu++)
    {
        //get type of orbital mu (s or p)
        type_mu = orbital_type[mu];

        //get valence charge of nucleus
        if (z_vals(atom_map[mu]) != 1)
        {
            z_a = z_vals(atom_map[mu]) - 2;
        }

        else 
        {
            z_a = 1;
        }
        

        for (int nu = 0; nu < n_basis; nu++)
        {
            //diagonal elements
            if (mu == nu)
            {
                //get paramters for orbital mu, s or p
                if (type_mu == 'S')
                {
                    constant = - s_parameters[orbital_species[mu]];
                }
                else if (type_mu == 'P')
                {
                    constant = - p_parameters[orbital_species[mu]];
                }

                //calculate sum in the expression for diagonal fock matrix elements
                sum=0;
                for (int atom_b = 0; atom_b < n_atoms; atom_b++)
                {
                    //get valence charge of nucleus
                    if (z_vals(atom_b) != 1)
                    {
                        z_b = z_vals(atom_b) - 2;
                    }

                    else 
                    {
                        z_b = 1;
                    
                    }

                    if (atom_map[mu] != atom_b)
                    {
                        sum += (z_b) * gamma_matrix(atom_map[mu], atom_b);
                    }
                }
                

                matrix(mu, mu) = constant - (z_a - (1.0/2.0)) * gamma_matrix(atom_map[mu], atom_map[mu]) - sum;
                                            
            }

            else 
            {
                matrix(mu, nu) = (1.0/2.0) * (beta[orbital_species[mu]] + beta[orbital_species[nu]])*overlap_matrix(mu, nu);
                                
            }
        }
    }

    return matrix;
}


//calculate the coefficient matrix for molecular orbitals using the given fock matrix
arma::mat molecule::calculate_coefficient_matrix(arma::mat fock_matrix)
{
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, fock_matrix);

    return eigenvectors;
}

//calculate the density matrix
arma::mat molecule::calculate_density_matrix(arma::mat coeff_matrix, int occ_mos)
{
    return coeff_matrix.cols(0, occ_mos - 1) * coeff_matrix.cols(0, occ_mos - 1).t();
}

arma::mat molecule::calculate_total_density_matrix(arma::mat density_alpha, arma::mat density_beta, arma::vec & atomic_density)
{
    arma::mat density_total = density_alpha + density_beta;

    for (int atom = 0; atom < n_atoms; atom++)
    {
        atomic_density(atom) = 0;
    }
    
    for (int mu = 0; mu < n_basis; mu++)
    {
        atomic_density(atom_map[mu]) += density_matrix_total(mu, mu);
    }
    return density_total;

}


void molecule::scf_steps(bool verbose)
{
    //p_electrons = 1;
    //q_electrons = 1;

    double tolerance = 1e-6;
    int iter = 0;

    arma::mat p_alpha_old;
    p_alpha_old.resize(n_basis, n_basis);
    arma::mat p_beta_old;
    p_beta_old.resize(n_basis, n_basis);

    coefficient_matrix_alpha = calculate_coefficient_matrix(fock_matrix_alpha);
    coefficient_matrix_beta = calculate_coefficient_matrix(fock_matrix_beta);

    if (verbose)
    {
        std::cout << "Performing SCF Iterations" << std::endl;
    }

    do {
        if (verbose)
        {
            std::cout << "Iteration: " << iter << std::endl;
        }
        arma::vec eigenvalues;
        arma::mat eigenvectors;



        p_alpha_old = density_matrix_alpha;
        density_matrix_alpha = calculate_density_matrix(coefficient_matrix_alpha, p_electrons);

        p_beta_old = density_matrix_beta;

        if (q_electrons > 0)
        {
            density_matrix_beta = calculate_density_matrix(coefficient_matrix_beta, q_electrons);
        }

        density_matrix_total = calculate_total_density_matrix(density_matrix_alpha, density_matrix_beta, atomic_density);


        fock_matrix_alpha = calculate_fock_matrix(density_matrix_alpha);
        arma::eig_sym(eigenvalues, eigenvectors, fock_matrix_alpha);
        coefficient_matrix_alpha = eigenvectors;

        
        if (q_electrons > 0)
        {
            

            fock_matrix_beta = calculate_fock_matrix(density_matrix_beta);
            arma::eig_sym(eigenvalues, eigenvectors, fock_matrix_beta);
            coefficient_matrix_beta = eigenvectors;
        }
        
        else 
        {
            density_matrix_total = density_matrix_alpha;

        }

        iter++;
       

    } while ((arma::norm(density_matrix_alpha - p_alpha_old, 2) > tolerance) || (arma::norm(density_matrix_beta - p_beta_old, 2) > tolerance));

}

double molecule::calculate_energy(bool verbose)
{
    double sum_alpha = 0;
    double sum_beta = 0;
    double sum_AB = 0;
    double result = 0;
    double conversion = 27.211;

    //int p_electrons = 1;
    //int q_electrons = 1;

    //density_matrix_alpha = calculate_density_matrix(coefficient_matrix_alpha, p_electrons);
    //density_matrix_beta = calculate_density_matrix(coefficient_matrix_beta, q_electrons);


    for (int mu = 0; mu < n_basis; mu++)
    {
        for (int nu = 0; nu < n_basis; nu++)
        {
            sum_alpha += density_matrix_alpha(mu, nu) * (core_hamiltonian(mu, nu) + fock_matrix_alpha(mu, nu));
            sum_beta +=  density_matrix_beta(mu, nu) * (core_hamiltonian(mu, nu) + fock_matrix_beta(mu, nu));

        }
    }

    double z_a;
    double z_b;
    for (int atom_a = 0; atom_a < n_atoms; atom_a++)
    {
        for (int atom_b = atom_a + 1; atom_b < n_atoms; atom_b++)
        {   
            if (z_vals(atom_a) != 1)
            {
                z_a = z_vals(atom_a) - 2;
            }

            else {
                z_a = 1;
            }
            
            if (z_vals(atom_b) != 1)
            {
                z_b = z_vals(atom_b) - 2;
            }

            else {
                z_b = 1;
            }
           
            sum_AB += (z_a * z_b) / arma::norm(coords.row(atom_a) - coords.row(atom_b), 2);
        }
    }

    result = ((1.0/2.0) * sum_alpha + (1.0/2.0) * sum_beta + sum_AB*conversion);
    if (verbose)
    {
        std::cout << "Nuclear Repulsion Energy: " << sum_AB*conversion << std::endl;
        std::cout << "Electron Energy: " << ((1.0/2.0) * sum_alpha + (1.0/2.0) * sum_beta)<< std::endl;
        std::cout << "Total Energy: " << result << std::endl;
    }
    

    return result;
}

void molecule::print_info()
{
    std::cout << "Gamma Matrix" << std::endl;
    std::cout << gamma_matrix << std::endl;

    std::cout << "Overlap Matrix" << std::endl;
    std::cout << overlap_matrix << std::endl;

    std::cout << "Coefficient Matrix Alpha" << std::endl;
    std::cout << coefficient_matrix_alpha << std::endl;

    std::cout << "Coefficient Matrix Beta" << std::endl;
    std::cout << coefficient_matrix_beta << std::endl;

    std::cout << "Density Matrix Alpha" << std::endl;
    std::cout << density_matrix_alpha << std::endl;

    std::cout << "Density Matrix Beta" << std::endl;
    std::cout << density_matrix_beta << std::endl;

    std::cout << "Atomic Density Vector" << std::endl;
    std::cout << atomic_density << std::endl;

    std::cout << "Fock Matrix Alpha" << std::endl;
    std::cout << fock_matrix_alpha << std::endl;

    std::cout << "Fock Matrix Beta" << std::endl;
    std::cout << fock_matrix_beta << std::endl;

    std::cout << "Core Hamiltonian" << std::endl;
    std::cout << core_hamiltonian << std::endl;
    
    std::cout << "Energy" << std::endl;
    std::cout << calculate_energy(true) << std::endl;
    std::cout << "------------------------------" << std::endl;
}



