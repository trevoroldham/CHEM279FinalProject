### Trevor Oldham
### CHEM 279
### UC Berkeley CoC

# Problem Set 5: Evaluate the analytic gradient of SCF energy

The code in this repository contains the functions and classes to perform the SCF iterations on a molecule using the CNDO/2 method. The analytical gradient is calculated using the equations from the problem statement. Input molecules are provided in /Molecules which specify the coordinates of the atoms, and the mutliplicity must be known by the user in order to provide both *p* and *q* to the constructor.

The project structure contains the following files and folders:

    /Molecules  -contains the molecules in input form
    atomic_orbital.cpp  -contains the atomic orbital class logic
    integrator.cpp  -contains the integration tools used in util functions
    main.cpp -contains the files for running test cases and optimizing bond length
    molecule.cpp    -contains the molecule class definition, functions and attributes
    primitive_guassian.cpp  -contains the primitive guassian class definition
    util.cpp    -contains mathematical functions used in molecule and integrator
    gradient.cpp   - contains the functions and logic used in evaluating the gradient


Using a closed or open shell molecule, the program can read the coordinates from a file in /Molecules and then perform the SCF minimization before printing relevant parameters as output including the density matrices, fock matrices, hamiltonian, and ground state energy. The functions in *gradient.cpp* are used to evaluate the gradient of a molecule after minimization. Test cases can be found in /Tests.

To compile use the command below:

    g++ -o main -std=c++17 main.cpp util.cpp integrator.cpp atomic_orbital.cpp primitive_guassian.cpp molecule.cpp gradient.cpp -DARMA_DONT_USE_WRAPPER -I /Users/trevor/Downloads/armadillo-12.6.2/include -framework Accelerate


