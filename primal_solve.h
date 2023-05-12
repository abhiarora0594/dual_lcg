#ifndef PRIMAL_SOLVE_H
#define PRIMAL_SOLVE_H

#include "dual_solve.h"

class Dual_Solve;
class Fem;
class Grid;
class Common_Utilities;

class Primal_Solve
{
private:

    PetscInt Istart,Iend;

public:
	Primal_Solve(Fem *fem, Grid *grid, Dual_Solve *dual_solve);

	~Primal_Solve();

    Fem *fem;
    Grid *grid;
    Dual_Solve *dual_solve;   
    Common_Utilities *common_utilities;

    int dof_per_node; // dof per node

    int dof_per_el; // dof per element

    double lambda_1; // stretch imposed along director

	double lambda_2; // stretch imposed in perp direction to director

    int count;

    int step_no;

    double residual_conv;

    double norm_rhs;

    int tdof;

    int tdof_fdef;

    MyMat <double> beta_c; // constraints on primal pde's

    MyMat <double> x_def; // x_def (deformed configuration)

    MyMat <double> F_def; // F_def (deformation gradient tensor)

    // mat and vecs for global newton solve (primal)
    Vec rhs_vec;
    Vec primal_vec;
    Vec delta_primal_vec;
    Mat K_primal;

    // Ksp solver, preconditioners for global newton solve (primal)
	KSP ksp_primal;
	PC pc_primal;
	Mat FF_primal; 

    // scattering for copying primal vector
    VecScatter 	ctx_primal;
	Vec  primal_vec_SEQ;

    // mat and vecs for l2 projection of F = grad x
    Vec fdef_vec;
    Vec F_fdef;
    Mat K_fdef;

    // Ksp solver, preconditioners for l2 projection of F = grad x
	KSP ksp_fdef;
	PC pc_fdef;
	Mat FF_fdef; 

    // scattering for copying f_def
    VecScatter 	ctx_fdef;
	Vec  fdef_vec_SEQ;

    // useful variables
    double *primal_sol;
    int *ind_set_primal;

    double *f_def_sol;
    int *ind_set_fdef;

    // methods
    PetscErrorCode set_coupling(Common_Utilities *);

    PetscErrorCode global_newton_raphson_primal(); 

	PetscErrorCode global_newton_solve();

    PetscErrorCode global_newton_update();

    PetscErrorCode elem_stiff_force_vector_global_newton(MyMat<double> , MyMat<double> ,
                                                            MyMat<double> , MyVec<double> ,
                                                            MyMat <double> ,
                                                            double * , double ** );

    PetscErrorCode copying_petsc_primal_vector_to_stdvector();

    PetscErrorCode initial_condition_for_primal_solve();

    PetscErrorCode deformation_gradient_l2_projection();

    PetscErrorCode elem_deformation_gradient_l2_projection_stiff_force(MyMat <double> ,
                                                MyVec <double> , double **, double *);

};


#endif