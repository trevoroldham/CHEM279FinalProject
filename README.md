### Trevor Oldham
### CHEM 279
### UC Berkeley CoC

# Final Project: CNDO/2 Method with DIIS Implementation and BFGS Optimization

The code in this repository contains the functions and classes to perform the SCF iterations on a molecule using the CNDO/2 method. The analytical gradient is calculated using the equations from the problem statement. Input molecules are provided in /Molecules which specify the coordinates of the atoms, and the mutliplicity must be known by the user in order to provide both *p* and *q* to the constructor.

The project structure contains the following files and folders:

    /Molecules  -contains the molecules in input form
    /Tests      -contains the test case degugging output for two sample molecules
    atomic_orbital.cpp  -contains the atomic orbital class logic
    integrator.cpp  -contains the integration tools used in util functions
    main.cpp -contains the files for running test cases and optimizing bond length using FP method and BFGS
    molecule.cpp    -contains the molecule class definition, functions and attributes
    primitive_guassian.cpp  -contains the primitive guassian class definition
    util.cpp    -contains mathematical functions used in molecule and integrator
    gradient.cpp   - contains the functions and logic used in evaluating the gradient
    SCF.cpp -contains the three algorithms for performing self consistent field minimization
    BFGS.cpp -contains the BFGS optimization
    Presentation.pdf - document with the presentation slides


Using a closed or open shell molecule, the program can read the coordinates from a file in /Molecules and then perform the SCF minimization before printing relevant parameters as output including the density matrices, fock matrices, hamiltonian, and ground state energy. The functions in *gradient.cpp* are used to evaluate the gradient of a molecule after minimization. Test cases can be found in /Tests.

To compile use the command below:

    g++ -o main -std=c++17 main.cpp util.cpp integrator.cpp atomic_orbital.cpp primitive_guassian.cpp molecule.cpp gradient.cpp SCF.cpp BFGS.cpp -DARMA_DONT_USE_WRAPPER -I /Users/trevor/Downloads/armadillo-12.6.2/include -framework Accelerate

### Direct Inversion of the Iterative Subspace

The *SCF* namespace contains three algorithms used to solve for the molecular coefficients in the CNDO/2 method. The first algorithm is a fixed point iteration that we have demonstrated in class. The second method is a mixture between fixed point iteration and DIIS, where the iterations continue until there are *m* entries in the fock matrix history, after which the algorithm switches to the DIIS method to extrapolate the new Fock matrix from a linear combination of the historical matrices. The final implementation is a pure DIIS method which begins solving the DIIS equations on the first iteration. The user can specify the desired tolerance and the number of iterations to save in the history.

More information about the algorithm can be found at the following reference:

https://en.wikipedia.org/wiki/DIIS

### Broyden–Fletcher–Goldfarb–Shanno Algorithm

The BFGS method is used to optimize the geometry of the molecule object. This method the gradient computed with the code in the *gradient* namespace and then a line search is performed to find the optimal step in direction of the gradient that falls under the Wolfe conditions. A vector *s* and a vector *y* are computed and then the Hessian is updated using the previous value of the Hessian as well as the values of *s* and *y*. This algorithm is still in the debugging process. The *main* file contains a brief example of the BFGS algorithm performed on the C2H4 molecule with verbose output detailing each step in the process. The algorithm does not converge until the maximum number of iterations is reached.

More information about the algorithm can be found at the following reference:

https://en.wikipedia.org/wiki/Broyden-Fletcher-Goldfarb-Shanno_algorithm


 

