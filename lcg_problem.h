#ifndef LCG_PROBLEM_H
#define LCG_PROBLEM_H

// list of constants
#define PI 3.141592653589793
#define LAMBDA 1.2
#define NU 1.5
#define T_START 0.0
#define T_FINAL 1.0
#define TOTAL_STEPS 5
#define STEP_NO 0
#define TOL_GLOBAL 1e-10
#define TOL_LOCAL 1e-6
#define PX 1e10
#define PF 1e10
#define Nprint 1
#define EPS 1.0
#define A_CONST 1e-4
#define PEN 1
#define EPS1 0.0
#define EPS2 0.0
#define EPS3 1.0
#define EPS3_res 1.0
#define PP 2
#define ONE_BY_EPS 0.0
#define PEN_CURL 0.0

enum MESH_INPUT
{
	cylindrical = 0,
	spherical,
	trelis_mesh,
	plane
};

enum PROBLEM_TYPE
{
	inverse = 0,
	forward
};

enum RUN_TYPE
{
	initial = 0,
	restart,
	restart_gp
};

enum DUAL_IC
{
	zero = 0,
	non_zero
};

enum ALGORITHM
{
	gradient_flow = 0,
	newton_raphson
};

MESH_INPUT mesh_input = trelis_mesh;
PROBLEM_TYPE problem_type = inverse;
RUN_TYPE run_type = initial;
ALGORITHM algorithm = gradient_flow;
DUAL_IC dual_ic = zero;


template <typename T>
using MyMat = std::vector<std::vector<T>>;

template <typename T>
using MyVec = std::vector<T>;


// Global variables
PetscInt 		nx=50,ny=50; // elements along x and y

PetscReal 		l1=1.0,l2=1.0; // size of the domain

PetscMPIInt    	rank; /* processor rank */

PetscMPIInt    	size; /* size of communicator */

// matrices and vector for solving curl part of A
Mat K_xi;
Vec F_xi;
Vec xi_vec;

// Ksp solver, preconditioners
KSP ksp_xi;
PC pc_xi;
Mat FF_xi;

// scattering data context to sequential vector
VecScatter 	ctx_xi;
Vec  xi_vec_SEQ;

// vectors for eigenvalue solvers
// Vec eigvAr, eigvAi;

// EPS eps; // eigenvalue solver context
// RG	rg; // region object

// PetscScalar kr, ki;
// PetscReal	error;

//signum function
inline int sgn(double val){
	return (0.0 < val) - (val < 0.0);
}

inline int sgn(int val){
	return (0 < val) - (val < 0);
}

#endif