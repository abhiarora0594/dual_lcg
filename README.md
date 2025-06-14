# dual_lcg
This file contains the code for the implementation of the dual variational principle
for inverse and forward design problems in the actuation of soft membranes.

# libraries
PETSc (MPI)

# code structure 
1. lcg_problem is the main file
2. Parameters to run in lcg_problem.h file
3. Grid class handles mesh
4. FEM class handles shape functions and its gradients
5. Dual_Solve is the main class
6. Common_Utilities class has some utility functions


# before running
1. create outputs directory --> stores files for outputting
2. create a restart directory --> stores files for restarting the code
