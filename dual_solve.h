#ifndef DUAL_SOLVE_H
#define DUAL_SOLVE_H


class Dual_Solve
{

private:
	double lambda_1; // stretch imposed along director

	double lambda_2; // stretch imposed in perp direction to director

	PetscInt Istart,Iend;

public:
	Dual_Solve(Fem *fem, Grid *grid, double lambda, double nu, \
					double t_start, double t_final, int total_steps, \
					int step_no, double tol_global, double tol_local, \
					double px, double pf);

	~Dual_Solve();

	Fem *fem; // pointer to Fem class

	Grid *grid; // pointer to Grid class

	int dim; // inverse or forward

	int F_DOF;

	double dt_step;

	double t_start;

	double t_step;

	double t_final;

	int step_no;

	int total_steps;

	int count;

	double residual_conv;

	double tol_global;

	double tol_local;

	double px,pf;

	Vec dual_vec;

	Vec delta_dual_vec;

	Vec rhs_vec;

	Mat K_dual;

	Vec rhs_vec_final;

	Vec mass_vec;

	VecScatter 	ctx_dual;

	Vec  dual_vec_SEQ;

	MyMat<double> II;

	MyMat<double> Id;

	MyMat <double> beta; // beta_1 and beta_2 (cartesian)

	MyMat <double> A_dual; // A dual field (c^k and E_gamma basis)

	MyMat <MyVec<double>> x_gp; // deformed coordinates at gpts (cartesian)

	MyMat <MyVec<double>> F_gp; // F tensor in (c_k and E^gamma basis) at gpts

	MyVec <double>	norm_global_x; // global vector for norm_x

	MyVec <double>  norm_global_F; // global vector for norm_F

	MyVec <double> dual_S_mu1;

	MyVec <double> dual_S_mu2;

	MyVec <double> dual_S_F_dx;

	MyVec <double> dual_S_quad_H;

	MyVec <double> dual_S_norm;

	MyVec <double> qfactor;

	MyVec <double> compatibility_check;

	double dual_S_total;

	double qf_global;

	double d_dual_dot;

	double undef_global; // undef global area of memberane

	// initial condition fields 

	double **x_def; // deformed co-ordinates

	MyMat <MyVec<double>> x_0_gp; // base state deformed co-ordinates

	MyMat <MyVec<double>> F_0_dual_gp; // F tensor for base state in (c^k and E_gamma basis) at gpts

	double *z_def; // z deformed coordinates.

	int tdof_z;

	int *ind_set_z; // dof nos for z deformed coordinates

	double *xi_dual; // xi- dual field

	int tdof_xi;

	int *ind_set_xi; // dof nos for xi_curl_A vector

	double *x_def_vec;

	double norm_rhs;
	// l2 projection for deformed shape

	// matrices and vector for solving l2 projection of deformed shape
	Mat K_xdef;

	Vec F_xdef;

	Vec xdef_vec;

	// matrices and vector for local newton solve
	Mat Ke_local;

	Vec Fe_local;

	Vec local_vec;

	// Ksp solver, preconditioners for deformed shape projection
	KSP ksp_xdef;

	PC pc_xdef;

	Mat FF_xdef;

	// Ksp solver, preconditioners for global newton solve
	KSP ksp_dual;

	PC pc_dual;

	Mat FF_dual;

	// Ksp solver, preconditioners for local newton solve
	KSP ksp_local;

	PC pc_local;

	Mat FF_local;

	// scattering data context to sequential vector
	VecScatter 	ctx_xdef;
	Vec  xdef_vec_SEQ;

	// gradient flow methods
	PetscErrorCode force_vector_assembly();

	PetscErrorCode gradient_flow();

	PetscErrorCode elem_force_vector(MyMat <double> , MyMat <double> , 
										 MyMat <double> , MyVec <double> , double*, 
										 MyMat <double> ,  MyMat <double> );

	PetscErrorCode elem_mass_vector(MyVec <double>, double*);

	PetscErrorCode picards_iteration_for_F(int , MyMat <double> , MyMat <double> ,
											MyMat<double>,	MyMat <double>, MyMat <double>, MyMat <double> &);

	PetscErrorCode local_newton_for_F(int , MyMat <double> , MyMat <double> ,
										MyMat <double> , MyMat <double> ,
										MyMat <double> , MyMat <double> &, int *);

	PetscErrorCode get_primal_fields(MyMat <double> , MyMat <double> ,
										MyMat <double>, MyMat <double>, MyMat <double>,
										MyMat <double>, MyMat <double> &, 
										MyMat <double> , MyMat <double> &);


	PetscErrorCode explicit_update();

	PetscErrorCode boundary_conditions();

	PetscErrorCode get_dual_l2_norm();

	PetscErrorCode get_elem_dual_l2_norm(MyMat <double> , MyMat <double> , 
									 MyMat <double> ,  MyVec<double> , 
									 MyMat <double> , MyMat <double> ,
									 MyMat <double> , MyMat <double> ,
									 double *, double *, 
									 double *, double *,
									 double *, double *,
									 double *, double *);

	// newton methods
	PetscErrorCode global_newton_raphson(); 

	PetscErrorCode global_newton_solve();

	PetscErrorCode global_newton_update();

	PetscErrorCode global_newton_solve_bcs(Mat &, Vec &, Vec &);

	PetscErrorCode elem_stiff_force_vector_global_newton(MyMat <double> , MyMat <double> , 
											 MyMat <double> ,  MyVec<double> , double* , 
											 double **, MyMat <double> ,  MyMat <double> );

	// helper methods 
	PetscErrorCode get_derivative_of_mus(MyMat <double>, double, double,
											MyMat <double>,	MyMat <double>, 
											int, int, double *, double *);

	PetscErrorCode second_derivative_of_invariants(MyVec <double>, MyVec <double> ,
													double , double , double , MyMat <double>,
													MyMat <double> , MyMat <double> , MyMat <double> ,
													MyMat <double> &, MyMat <double> &);

	PetscErrorCode first_derivative_of_invariants(MyVec <double> , MyVec <double> ,
													MyMat <double> , MyMat<double> ,
													MyVec <double> &, MyVec <double> &);

	inline PetscErrorCode get_eigenvectors_of_C(MyMat <double> , double , double , MyMat <double>,
												MyVec <double> &, MyVec <double> &);

	inline PetscErrorCode tranform_C_on_both_legs(MyMat <double> ,MyMat <double> ,
														MyMat <double> , MyMat <double> &);

	inline PetscErrorCode get_G_dual(int , MyMat <double>, MyMat<double> &);

	inline PetscErrorCode get_H_term(MyMat <double> , MyMat <double> ,
										int , int, double *);

	inline PetscErrorCode get_A_field_at_gp(MyMat <double> ,
												int , int , double *);

	inline PetscErrorCode get_beta_field_at_gp(MyMat <double> , int , 
												double *, double *);

	inline PetscErrorCode get_invariants(double , double , double *, double *);

	inline PetscErrorCode trace_of_C(MyMat <double> , double*);

	inline PetscErrorCode determinant_of_C(MyMat <double> , MyMat <double> , double *);

	inline PetscErrorCode determinant_of_F(MyMat <double> , MyMat <double> , double *);

	inline PetscErrorCode tranform_F_to_dual(MyMat<double> , MyMat<double> , MyMat<double> &);

	inline PetscErrorCode tranform_F_from_dual(MyMat<double>, MyMat<double> , MyMat<double> &);

	inline PetscErrorCode get_C(MyMat<double> , MyMat <double> , MyMat <double> &);

	inline PetscErrorCode get_div_A_field_at_gp(MyMat<double> , int , double *);

	template<class T>
	inline PetscErrorCode get_max_val(MyVec<T> v, T *val);

	// initial condition methods

	PetscErrorCode initial_condition_dual();

	PetscErrorCode initial_condition_primal_to_dual();

	PetscErrorCode initial_guess_cylindrical();

	PetscErrorCode initial_guess_trelis();

	PetscErrorCode initial_guess_plane();

	PetscErrorCode initial_guess_hat_shape();

	PetscErrorCode get_x_F_at_gps_initial_guess();

	PetscErrorCode get_x_F_at_gps_restart_guess();

	PetscErrorCode solve_grad_part_A_dual();

	PetscErrorCode solve_curl_part_A_dual();

	PetscErrorCode solve_dual_fields_from_primal();

	PetscErrorCode elem_grad_A_stiff_force(MyMat<double> , MyVec<double> , double [][8], double *);

	PetscErrorCode grad_A_solve_bcs(Mat &, Vec &, Vec &);

	PetscErrorCode A_solve_bcs(Mat &, Vec &, Vec &);

	PetscErrorCode dual_solve_bcs(Mat &, Vec &, Vec &);

	PetscErrorCode elem_L2_A_stiff_force(MyMat <double> ,MyMat <double> , MyVec<double> ,
											MyMat <double> , MyMat<double> ,
											double ke[][16], double *fe);

	PetscErrorCode elem_L2_dual_stiff_force(MyMat <double> , MyMat <double> ,MyMat <double> , MyMat <double> , 
												MyVec<double> ,	MyMat <double> , double **, double *);

	PetscErrorCode compute_L2_norm_initial_guess(); 

	PetscErrorCode elem_compute_L2_norm_initial_guess(MyMat<double> , MyMat<double>, 
														MyMat<double> , MyMat<double> ,
														MyMat<double> , MyMat<double> , MyMat<double> ,
														MyVec<double> ,double *, double *);

	PetscErrorCode deformed_shape_l2_projection();

	PetscErrorCode elem_deformed_geometry_solve_stiff_force(MyMat<double> , MyMat<double>, MyVec<double> ,
															double **, double *);

	PetscErrorCode area_calculation();
	
	// petsc methods
	PetscErrorCode mat_create_petsc(Mat &, int, int, int);

	PetscErrorCode vec_create_petsc(Vec &, int);

	PetscErrorCode ksp_mumps_solver_petsc(KSP &, Mat &, PC &, Mat &);

	PetscErrorCode delete_petsc_objects();

	PetscErrorCode copying_petsc_dual_vector_to_stdvector();

	PetscErrorCode copying_dual_to_petsc_vector();

	// outputs
	PetscErrorCode vtk_write();

	PetscErrorCode restart_write();

	// small liner system solve
	PetscErrorCode mysolve_local_linear_system(int , MyMat <double> ,
												MyVec <double> , MyVec <double> &);

	PetscErrorCode petsc_local_linear_system(int , MyMat <double> ,
												MyVec <double> , MyVec <double> &);												

};


#endif