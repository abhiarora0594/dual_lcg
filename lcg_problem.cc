static char help[] = "Membrane problem of finding an imprinting pattern\n\
for a given reference membrane to achieve a desired target membrance under photo/thermal stimulus.\n\
Input parameters include:\n\
  -nx  <mesh_x>       	: number of elements in x-direction\n\
  -ny  <mesh_y>       	: number of elements in y-direction\n\
  -l1 <length_x>    : domain length in x-direction\n\
  -l2 <length_y>    : domain length in y-direction\n\n";

#include <petscksp.h>	 // Include "petscksp.h" so that we can use KSP solvers
// #include <slepceps.h>  // Include "slepsceps.h" to solve eigenvalue problem
#include <time.h>			 // Import this package to print time to stdout
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <cmath>  
#include "lcg_problem.h"  // Include parameters and macros
#include "grid.cc" // grid class for mesh 
#include "fem.cc" // fem class for shape fns, derivaties
#include "primal_solve.cc"
#include "common_utilities.cc"
#include "dual_solve.cc"  // main class for gradient flow or Newton raphson of dual variational principle


using namespace std;

int main(int argc,char **args)
{

	//petsc variables
	PetscErrorCode ierr;

	// initialize PetSc and take in input arguments
	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
	ierr = PetscOptionsGetReal(NULL,NULL,"-l1",&l1,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-l2",&l2,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&(::rank));CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&(::size));CHKERRQ(ierr);

	time_t T=time(NULL);
	struct tm tm=*localtime(&T);

	// initiate grid class here 
	Grid *grid = new Grid(nx,ny,l1,l2);

	// creating mesh
	ierr = grid->mesh_data_read(); CHKERRQ(ierr);

	// get minimum h from mesh
	ierr = grid->get_h_min(); CHKERRQ(ierr);

	// initiate fem class here
	Fem *fem = new Fem(grid,4,4);

	//shape funcs and its derivatives
	ierr = fem->jacobian_convected_basis(); CHKERRQ(ierr);

	Dual_Solve *dual_solve = new Dual_Solve(fem,grid,LAMBDA,NU,T_START,T_FINAL,
												TOTAL_STEPS,STEP_NO,TOL_GLOBAL,TOL_LOCAL,PX,PF);

	std::cout << "All classes constructed" << std::endl;

	// contruct the primal solve class
	ierr = dual_solve->contruct_primal_solve_class(); CHKERRQ(ierr);

	// construct common utilities class
	ierr = dual_solve->contruct_common_utilities_class(); CHKERRQ(ierr);

	// set coupling
	ierr = dual_solve->set_coupling(); CHKERRQ(ierr);

	// dual solve
	if (algorithm == gradient_flow){
		ierr = dual_solve->gradient_flow(); CHKERRQ(ierr);	
	}
	else if (algorithm == newton_raphson){
		ierr = dual_solve->global_newton_raphson(); CHKERRQ(ierr);
	}
	

	// delete objects
	delete dual_solve;
	delete fem;
	delete grid;

	std::cout << "All classes destructed" << std::endl;
	
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = PetscFinalize(); CHKERRQ(ierr);
	
	return ierr;

}