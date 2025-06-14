#include "dual_solve.h"
#include "common_utilities.h"
#include "primal_solve.h"


Dual_Solve::Dual_Solve(Fem *fem, Grid *grid, double lambda, double nu, \
						double t_start, double t_final, int total_steps, \
						int step_no, double tol_global, double tol_local, \
						double px, double pf)
{

	double large_root;
	double small_root;

	this->fem = fem;

	this->grid = grid;

	if (problem_type == inverse){
		this->dim = 2;
		large_root = 1.0/std::pow(lambda,-nu);
		small_root = 1.0/lambda;
		this->lambda_1 = 1.0 + (4.0+1.0)*(large_root-1.0)/5.0;
		this->lambda_2 = 1.0 - (4.0+1.0)*(1.0-small_root)/5.0;
	}
	else if (problem_type == forward){
		this->dim = 3;
		this->lambda_1 = std::sqrt(2.0); //lambda;
		this->lambda_2 = 1.0; // std::pow(lambda,-nu);
	}

	this->F_DOF = 2*this->dim;

	this->Istart = (::rank)*(grid->nel/(::size)) + ((grid->nel%(::size)) < (::rank) ? (grid->nel%(::size)) : (::rank));

  	this->Iend = Istart + grid->nel/(::size) + ((grid->nel%(::size)) > (::rank));

  	std::cout << " lambda_1 is " << this->lambda_1 << " and lambda_2 is " << this->lambda_2 << std::endl;

  	std::cout << " dim is " << this->dim << std::endl;

  	this->t_start = t_start;

  	this->t_final = t_final;

  	this->total_steps = total_steps;

  	this->t_step = this->t_start;

  	this->tol_global = tol_global;

  	this->tol_local = tol_local;

  	this->step_no = step_no;

  	this->px = px;

  	this->pf = pf;

  	this->count = 0;

  	II = MyMat<double> (dim,MyVec<double>(dim));

	Id = MyMat<double> (2,MyVec<double>(2));

	beta = MyMat <double> (grid->ndof,MyVec<double>(2)); 
	// beta_1 and beta_2 (cartesian)

	A_dual = MyMat <double> (grid->ndof,MyVec<double>(grid->A_DOF)); 
	// A dual field (c^k and E_gamma basis)

	beta_rhs = MyMat <double> (grid->ndof,MyVec<double>(2)); 
	// beta_1 and beta_2 (cartesian)

	A_dual_rhs = MyMat <double> (grid->ndof,MyVec<double>(grid->A_DOF)); 
	// A dual field (c^k and E_gamma basis)

	x_gp = MyMat <MyVec<double>> (grid->nel,MyMat<double>(fem->ngp,MyVec<double>(dim))); 
	// deformed coordinates at gpts (cartesian)

	F_gp = MyMat <MyVec<double>> (grid->nel,MyMat<double>(fem->ngp,MyVec<double>(F_DOF))); 
	// F tensor in (c_k and E^gamma basis) at gpts	

	x_0_def = MyMat <double> (grid->ndof,MyVec<double>(this->dim));
	// x deformed coordinates

	// base state deformed coordinates at gpts
	x_0_gp = MyMat <MyVec<double>> (grid->nel,MyMat<double>(fem->ngp,MyVec<double>(dim)));

	// base state F tensor in dual basis at gpts
	F_0_dual_gp = MyMat <MyVec<double>> (grid->nel,MyMat<double>(fem->ngp,MyVec<double>(F_DOF))); 

	// base state F tensor in real basis at gpts
	F_0_gp = MyMat <MyVec<double>> (grid->nel,MyMat<double>(fem->ngp,MyVec<double>(F_DOF))); 

	Eig_1_gp = MyMat <MyVec <double>> (grid->nel, MyMat <double>(fem->ngp,MyVec <double>(grid->X_DOF)));

	Eig_2_gp = MyMat <MyVec <double>> (grid->nel, MyMat <double>(fem->ngp,MyVec <double>(grid->X_DOF)));

	x_def = new double* [grid->ndof];
	for (int i=0;i<grid->ndof;i++)
		x_def[i]= new double [this->dim];

	x_def_vec = new double [grid->ndof*this->dim];

	// L2 norm x
	norm_global_x = MyVec <double> (grid->nel);

	// L2 norm F
	norm_global_F = MyVec <double> (grid->nel);

	dual_S_mu1 = MyVec <double> (grid->nel);

	dual_S_mu2 = MyVec <double> (grid->nel);

	dual_S_F_dx = MyVec <double> (grid->nel);

	dual_S_quad_H = MyVec <double> (grid->nel);

	dual_S_norm = MyVec <double> (grid->nel);

	qfactor = MyVec <double> (grid->nel);

	compatibility_check = MyVec <double> (grid->nel);

	// z- deformed co-ordinates
  	this->tdof_z = this->dim*grid->ndof;

  	z_def = new double [this->tdof_z];

  	ind_set_z = new int [this->tdof_z];

  	// xi dual vector 
  	this->tdof_xi = grid->dof_per_node*grid->ndof;

  	xi_dual = new double [this->tdof_xi];

  	ind_set_xi = new int [this->tdof_xi];


  	for (int i=0;i<this->dim;i++){
  		for (int j=0;j<this->dim;j++){

  			if (i==j)
  				this->II[i][j] = 1.0;
  			else
  				this->II[i][j] = 0.0;
  		}
  	}

  	for (int i=0;i<2;i++){
  		for (int j=0;j<2;j++){

  			if (i==j)
  				this->Id[i][j] = 1.0;
  			else
  				this->Id[i][j] = 0.0;
  		}
  	}


  	// index set for z (grad part of A_dual)
  	for (int i=0;i<grid->ndof;i++){
  		for (int d=0;d<this->dim;d++){

  			ind_set_z[this->dim*i+d] = this->dim*i+d;
  		}
  	}

  	// index set for xi (curl part of A_dual)
  	for (int i=0;i<grid->ndof;i++){
  		for (int d=0;d<grid->dof_per_node;d++){

  			ind_set_xi[grid->dof_per_node*i+d] = grid->dof_per_node*i+d;
  		}
	}

}

Dual_Solve::~Dual_Solve()
{

	// deleting variables and classes
	// petsc objects deleted separetely in a function

	delete[] ind_set_z;

	delete[] z_def;

	delete[] ind_set_xi;

	delete[] xi_dual;

	for (int i=0; i<grid->ndof; i++)
		delete[] x_def[i];

	delete[] x_def;

	delete[] x_def_vec;

	// deleting classes
	delete primal_solve;
	delete common_utilities;
}

PetscErrorCode Dual_Solve::contruct_primal_solve_class()
{
	// constructing primal solve class
	this->primal_solve = new Primal_Solve(this->fem,this->grid,this);

	return 0;
}


PetscErrorCode Dual_Solve::contruct_common_utilities_class()
{
	// constructing common_utilities class
	this->common_utilities = new Common_Utilities(this->fem,this->grid,this);

	return 0;
}

PetscErrorCode Dual_Solve::set_coupling()
{
	// primal solve couple
	primal_solve->set_coupling(common_utilities);

	return 0;
}

PetscErrorCode Dual_Solve::gradient_flow()
{

	PetscErrorCode ierr;

	double qf_max = 0.0;

	std::ofstream myfile;

	if (run_type == initial){
		myfile.open("./outputs/timepoints.txt", std::ios::out);
	}
	else if (run_type == restart || run_type == restart_gp){
		myfile.open("./outputs/timepoints.txt", std::ios::app);		
	}

	if (run_type == restart || run_type == restart_gp){

		std::ifstream time_stamp;

		time_stamp.open("restart_t.txt", std::ios::in);

		time_stamp >> step_no >> t_step;

		time_stamp.close();

		t_start = t_step;

		std::cout << "step_no and t_step is " << step_no << " " << t_step << std::endl;

	}

	if (::rank == 0){

		myfile << "%%step_no" << " " << "dt_step" << " " << "t_step" << " " << 
					"dual_S_total"  << " " << "norm_rhs" << " " << "qf_max" << 
					" " << "qf_global" << " " << "dual_diff" << std::endl;
	}

	ierr = common_utilities->vec_create_petsc(rhs_vec,grid->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(rhs_vec_final,grid->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(dual_vec,grid->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(mass_vec,grid->tdof); CHKERRQ(ierr);

	// temporary variables
	int xdef_dof = grid->ndof*this->dim;

	std::cout << "xdef_dof is " << xdef_dof << std::endl;

	// ierr = common_utilities->mat_create_petsc(K_xdef,xdef_dof,30,25); CHKERRQ(ierr);
	ierr = common_utilities->mat_create_petsc(K_xdef,xdef_dof,40,35); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(F_xdef,xdef_dof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(xdef_vec,xdef_dof); CHKERRQ(ierr);


	try	
	{

		ierr = this->area_calculation(); CHKERRQ(ierr);

		std::cout<< "area undef is " << undef_global << std::endl;

		ierr = this->initial_condition_before_primal_solve(); CHKERRQ(ierr);

		if (run_type == initial){
			// ierr = primal_solve->global_newton_raphson_primal(); CHKERRQ(ierr);
		}

		ierr = this->initial_condition_after_primal_solve(); CHKERRQ(ierr);

		ierr = this->initial_condition_dual(); CHKERRQ(ierr);
		
		// ierr = this->initial_condition_primal_to_dual(); CHKERRQ(ierr);

		ierr = this->get_dual_l2_norm(); CHKERRQ(ierr);

		ierr = common_utilities->get_max_val(qfactor, &qf_max); CHKERRQ(ierr);

		residual_conv = 1.0;
	
		while ((t_step < t_final || residual_conv > tol_global) && step_no < total_steps)
		{

			if ( count == 0  || step_no == 10 || step_no%Nprint==0){
				ierr = this->deformed_shape_l2_projection(); CHKERRQ(ierr);
			}

			if ( count == 0 && ::rank == 0){
				ierr = common_utilities->vtk_write(); CHKERRQ(ierr);
				ierr = common_utilities->restart_write(); CHKERRQ(ierr);
			}

			ierr = this->force_vector_assembly(); CHKERRQ(ierr);

			ierr = this->explicit_update(); CHKERRQ(ierr);

			ierr = this->get_dual_l2_norm(); CHKERRQ(ierr);

			qf_global = std::sqrt(qf_global/undef_global);

			ierr = common_utilities->get_max_val(qfactor, &qf_max); CHKERRQ(ierr);

			t_step = t_step + dt_step;

			step_no = step_no + 1;

			count = count + 1;

			residual_conv = norm_rhs;

			if (::rank == 0){

				std::cout << "**********************************" << std::endl;
				std::cout << "step_no is " << step_no << std::endl;
				std::cout << "dt_step is " << dt_step << std::fixed << std::setprecision(15) << std::endl;
				std::cout << "t_step is " << t_step << std::endl;
				std::cout << "dual_S_total is " << dual_S_total  << std::endl;
				std::cout << "norm_rhs is " << norm_rhs << std::endl;
				std::cout << "qf_max is " <<  qf_max << std::endl;
				std::cout << "qf_global is " << qf_global << std::endl; 
				std::cout << "d_dual_dot is " << d_dual_dot << std::endl;
				std::cout << "**********************************" << std::endl;

				myfile << step_no << " " <<
							dt_step << " " <<  std::fixed << std::setprecision(15) <<
							t_step << " " << 
							dual_S_total  << " " << 
							norm_rhs << " " << 
							qf_max << " " << 
							qf_global << " " <<
							d_dual_dot << std::endl; 	
			}
			
			if ((step_no == 10 || step_no%Nprint==0) && ::rank==0){
				ierr = common_utilities->vtk_write(); CHKERRQ(ierr);
				ierr = common_utilities->restart_write(); CHKERRQ(ierr);	
			}
			

		} // while loop ends here

		if ((step_no%Nprint != 0) && ::rank==0){
			ierr = common_utilities->vtk_write(); CHKERRQ(ierr);
			ierr = common_utilities->restart_write(); CHKERRQ(ierr);	
		}

	}
	catch (const char *exception)
	{
		std::cerr << "Error: " << exception << "\n";
	}

	myfile.close();	

	return 0;
}


PetscErrorCode Dual_Solve::global_newton_raphson()
{

	PetscErrorCode ierr;

	double qf_max = 0.0;

	dt_step = 1.0;

	std::ofstream myfile;

	if (run_type == initial){
		myfile.open("./outputs/timepoints.txt", std::ios::out);
	}
	else if (run_type == restart || run_type == restart_gp){
		myfile.open("./outputs/timepoints.txt", std::ios::app);		
	}

	if (run_type == restart || run_type == restart_gp){

		std::ifstream time_stamp;

		time_stamp.open("restart_t.txt", std::ios::in);

		time_stamp >> step_no >> t_step;

		time_stamp.close();

		t_start = t_step;

		std::cout << "step_no is " << step_no << std::endl;

	}

	if (::rank == 0){

		myfile << "%%step_no" << " " << "dual_S_total"  << " " << "norm_rhs" 
				<< " " << "qf_max" << " " << "qf_global"  << std::endl;
	}	

	ierr = common_utilities->vec_create_petsc(rhs_vec,grid->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(dual_vec,grid->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(delta_dual_vec,grid->tdof); CHKERRQ(ierr);
	ierr = common_utilities->mat_create_petsc(K_dual,grid->tdof,150,150); CHKERRQ(ierr);
	ierr = VecSet(rhs_vec,0.0); CHKERRQ(ierr);
	ierr = VecSet(delta_dual_vec,0.0); CHKERRQ(ierr);

	// ierr = common_utilities->vec_create_petsc(local_vec,F_DOF); CHKERRQ(ierr);
	// ierr = common_utilities->vec_create_petsc(Fe_local,F_DOF); CHKERRQ(ierr);
	// ierr = common_utilities->mat_create_petsc(Ke_local,F_DOF,F_DOF,F_DOF); CHKERRQ(ierr);
	// ierr = VecSet(local_vec,0.0); CHKERRQ(ierr);
	// ierr = VecSet(Fe_local,0.0); CHKERRQ(ierr);

	// temporary variables
	int xdef_dof = grid->ndof*this->dim;

	std::cout << "xdef_dof is " << xdef_dof << std::endl;

	ierr = common_utilities->mat_create_petsc(K_xdef,xdef_dof,40,35); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(F_xdef,xdef_dof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(xdef_vec,xdef_dof); CHKERRQ(ierr);

	try	
	{

		ierr = this->area_calculation(); CHKERRQ(ierr);

		std::cout<< "area undef is " << undef_global << std::endl;

		ierr = this->initial_condition_before_primal_solve(); CHKERRQ(ierr); 

		if (run_type == initial){
			// ierr = primal_solve->global_newton_raphson_primal(); CHKERRQ(ierr);
		}

		ierr = this->initial_condition_after_primal_solve(); CHKERRQ(ierr); 

		ierr = this->initial_condition_dual(); CHKERRQ(ierr); 

		ierr = this->get_dual_l2_norm(); CHKERRQ(ierr);

		ierr = common_utilities->get_max_val(qfactor, &qf_max); CHKERRQ(ierr);

		residual_conv = 1.0;
	
		while (residual_conv > tol_global || step_no < TOTAL_STEPS)
		{

			if (count == 0){
				ierr = this->deformed_shape_l2_projection(); CHKERRQ(ierr);
			}

			// std::cout << " code is running fine before printing write files " << std::endl;

			if (count == 0 && ::rank == 0){

				ierr = common_utilities->vtk_write(); CHKERRQ(ierr);
				ierr = common_utilities->restart_write(); CHKERRQ(ierr);
			}

			// std::cout << " code is running fine after printing write files " << std::endl;

			ierr = this->global_newton_solve(); CHKERRQ(ierr);

			// std::cout << " code is running fine after global newton solve " << std::endl;

			ierr = this->global_newton_update(); CHKERRQ(ierr);

			ierr = VecSet(rhs_vec,0.0); CHKERRQ(ierr);
			
			ierr = VecSet(delta_dual_vec,0.0); CHKERRQ(ierr);

			ierr = MatSetOption(K_dual,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);

			ierr = MatZeroEntries(K_dual); CHKERRQ(ierr);

			ierr = this->get_dual_l2_norm(); CHKERRQ(ierr);

			qf_global = std::sqrt(qf_global/undef_global);

			ierr = common_utilities->get_max_val(qfactor, &qf_max); CHKERRQ(ierr);

			step_no = step_no + 1;

			count = count + 1;

			t_step = t_step + dt_step;

			residual_conv = norm_rhs;

			if (::rank == 0){

				std::cout << "**********************************" << std::endl;
				std::cout << "step_no is " << step_no << std::endl;
				std::cout << "dual_S_total is " << dual_S_total << std::fixed << std::setprecision(15) << std::endl;
				std::cout << "norm_rhs is " << norm_rhs << std::endl;
				std::cout << "qf_max is " <<  qf_max << std::endl;
				std::cout << "qf_global is " << qf_global << std::endl; 
				std::cout << "**********************************" << std::endl;

				myfile << step_no << " " <<
							dual_S_total  << " " << std::fixed << std::setprecision(15) <<
							norm_rhs << " " << 
							qf_max << " " << 
							qf_global << std::endl; 	
			}

			if (step_no == 10 || step_no%Nprint==0){
				ierr = this->deformed_shape_l2_projection(); CHKERRQ(ierr);
			}
			
			if ((step_no == 10 || step_no%Nprint==0) && ::rank==0){
				ierr = common_utilities->vtk_write(); CHKERRQ(ierr);
				ierr = common_utilities->restart_write(); CHKERRQ(ierr);	
			}
			

		} // while loop ends here

		if ((step_no%Nprint != 0) && ::rank==0){
			ierr = common_utilities->vtk_write(); CHKERRQ(ierr);
			ierr = common_utilities->restart_write(); CHKERRQ(ierr);	
		}

	}
	catch (const char *exception)
	{
		std::cerr << "Error: " << exception << "\n";
	}

	myfile.close();	

	return 0;
}

PetscErrorCode Dual_Solve::global_newton_update()
{
	PetscErrorCode ierr;

	double scal = 1.0;

	ierr = VecNorm(rhs_vec, NORM_INFINITY, &norm_rhs); CHKERRQ(ierr);

	ierr = VecAXPY(dual_vec,scal,delta_dual_vec); CHKERRQ(ierr);

	ierr = this->copying_petsc_dual_vector_to_stdvector(); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode Dual_Solve::explicit_update()
{	

	PetscErrorCode ierr;

	// double tmp_dt;

	double max_mass, min_mass;

	ierr = this->boundary_conditions(); CHKERRQ(ierr);

	ierr = VecPointwiseDivide(rhs_vec_final,rhs_vec,mass_vec); CHKERRQ(ierr);

	ierr = VecNorm(rhs_vec, NORM_INFINITY, &norm_rhs); CHKERRQ(ierr);

	ierr = VecNorm(rhs_vec_final, NORM_INFINITY, &d_dual_dot); CHKERRQ(ierr);

	ierr = VecNorm(mass_vec,NORM_INFINITY,&max_mass); CHKERRQ(ierr);

	// dt_step = std::min(0.1/norm_rhs, 10.0*grid->h_min);
	
	dt_step = 0.1/norm_rhs;
	// dt_step = 0.1;

	ierr = VecMin(mass_vec,NULL,&min_mass); CHKERRQ(ierr);

	// tmp_dt = dt_step; //from previous time step

	// // dt_step = std::max(100.0/norm_rhs,1.2*tmp_dt);

	// if (dual_diff < 1.1){
	// 	// tmp_dt = dt_step;
	// 	alpha = alpha + 5.0;
	// 	alpha = std::min(alpha,500.0);
	// 	dt_step = alpha/norm_rhs;
	// 	std::cout << "increasing alpha to " << alpha << " and dt_step = alpha/norm_rhs" << std::endl;
	// }
	// else {
	// 	alpha = 10.0;
	// 	dt_step = 10.0/norm_rhs;
	// 	std::cout << "cutting back to dt_step = 10.0/norm_rhs" << std::endl;
	// }

	// std::cout << "mass_vec max and min is " << max_mass << " " << min_mass << std::endl;

	ierr = VecAXPY(dual_vec,dt_step,rhs_vec_final); CHKERRQ(ierr);

	ierr = this->copying_petsc_dual_vector_to_stdvector(); CHKERRQ(ierr);

	return 0;
}


PetscErrorCode Dual_Solve::area_calculation()
{

	double tmp=0.0;

	for (int ie=0;ie<grid->nel;++ie){
		for (int i=0;i<fem->ngp;++i){

			tmp = tmp + fem->detJ_G[ie][i]*fem->wt[i];
		}
	}

	undef_global = tmp;

	return 0;
}


PetscErrorCode Dual_Solve::initial_condition_dual()
{	

	PetscErrorCode ierr;	

	int indc[grid->dof_per_node];
	double val[grid->dof_per_node];

	if (dual_ic == zero)
	{
		if (problem_type == inverse)
		{
			for (int i=0;i<grid->ndof; i++){

				beta[i][0] = 0.0;
				beta[i][1] = 0.0;

				A_dual[i][0] = 0.0;
				A_dual[i][1] = 0.0;
				A_dual[i][2] = 0.0;
				A_dual[i][3] = 0.0;
			}
		}
		else if (problem_type == forward)
		{
			for (int i=0;i<grid->ndof; i++){
				
				beta[i][0] = 0.0;
				beta[i][1] = 0.0;

				A_dual[i][0] = 0.0;
				A_dual[i][1] = 0.0;
				A_dual[i][2] = 0.0;
				A_dual[i][3] = 0.0;
				A_dual[i][4] = 0.0;
				A_dual[i][5] = 0.0;
			}
		}
	}
	else if (dual_ic == non_zero)
	{
		std::ifstream myfile;

		myfile.open("restart_dual.txt",std::ios::in);

		if (problem_type == inverse)
		{
			for(int i=0;i<grid->ndof;i++){
				myfile >> beta[i][0] >> beta[i][1] >>
							A_dual[i][0] >> A_dual[i][1] >> 
							A_dual[i][2] >> A_dual[i][3];	
			}
		}
		else if (problem_type == forward)
		{
			for(int i=0;i<grid->ndof;i++){
				myfile >> beta[i][0] >> beta[i][1] >>
							A_dual[i][0] >> A_dual[i][1] >> 
							A_dual[i][2] >> A_dual[i][3] >>
							A_dual[i][4] >> A_dual[i][5];	
			}
		}

		myfile.close();
	}


	for (int i=0;i<grid->ndof;i++){

		indc[0] = i*grid->dof_per_node;
		indc[1] = i*grid->dof_per_node + 1;

		val[0] = beta[i][0];
		val[1] = beta[i][1];

		for (int k1=0;k1<this->dim;k1++){
			for (int k2=0;k2<2;k2++){

				val[2+2*k1+k2] = A_dual[i][2*k1+k2];

				indc[2+2*k1+k2] = i*grid->dof_per_node + 2 + 2*k1 + k2;
			}
		}

		ierr = VecSetValues(dual_vec,grid->dof_per_node,indc,val,INSERT_VALUES); CHKERRQ(ierr);

		// final assembly of global vector (due to multiple proccess)
		ierr = VecAssemblyBegin(dual_vec); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(dual_vec); CHKERRQ(ierr);

	}

	return 0;
}

PetscErrorCode Dual_Solve::initial_condition_before_primal_solve()
{

	PetscErrorCode ierr;

	if (run_type == initial)
	{
		// deformed coordinates guess
		if (mesh_input == cylindrical){
			ierr = this->initial_guess_cylindrical(); CHKERRQ(ierr);
		}
		else if (mesh_input == plane){
			ierr = this->initial_guess_plane(); CHKERRQ(ierr);
		}
		else if (mesh_input == trelis_mesh || mesh_input == spherical){
			
			ierr = this->initial_guess_trelis(); CHKERRQ(ierr);

			// ierr = this->initial_guess_hat_shape(); CHKERRQ(ierr);
		}
		
		if (problem_type == inverse)
		{
			// base state same as initial state
			for(int i=0;i<grid->ndof;i++){

				x_0_def[i][0] = x_def[i][0];
				x_0_def[i][1] = x_def[i][1];	
			}
		}
		else if (problem_type == forward)
		{
			// base state same as initial state
			for(int i=0;i<grid->ndof;i++){

				x_0_def[i][0] = x_def[i][0];
				x_0_def[i][1] = x_def[i][1];
				x_0_def[i][2] = x_def[i][2];	
			}
		}
		
	}
	else if (run_type == restart)
	{

		std::ifstream myfile;

		myfile.open("restart_x_def.txt",std::ios::in);

		for(int i=0;i<grid->ndof;i++){

			myfile >> x_def[i][0] >> x_def[i][1];	
		}

		myfile.close();

	}

	return 0;
}

PetscErrorCode Dual_Solve::initial_condition_after_primal_solve()
{ 

	PetscErrorCode ierr;

	if (run_type == initial){
		// for (int i=0;i<grid->ndof;i++){

		// 	this->x_def[i][0] = primal_solve->x_def[i][0];
		// 	this->x_def[i][1] = primal_solve->x_def[i][1];
		// }
	}

	// get initial coordinates (x_gp) and deformation gradient tensor (F_gp) at gpts
	if (run_type == restart_gp)
	{
		ierr = this->get_x_F_at_gps_restart_guess(); CHKERRQ(ierr);
	}
	else if (run_type == initial || run_type == restart)
	{	
		ierr = this->get_x_F_at_gps_initial_guess(); CHKERRQ(ierr);	
	}

	return 0;
}


PetscErrorCode Dual_Solve::initial_condition_primal_to_dual()
{	

	PetscErrorCode ierr;

	// deformed coordinates guess
	if (mesh_input == cylindrical){
		ierr = this->initial_guess_cylindrical(); CHKERRQ(ierr);
	}
	else if (mesh_input == trelis_mesh || mesh_input == spherical){
		
		ierr = this->initial_guess_trelis(); CHKERRQ(ierr);
		// ierr = this->initial_guess_hat_shape(); CHKERRQ(ierr);
	}

	// get initial coordinates (x_gp) and deformation gradient tensor (F_gp) at gpts
	ierr = this->get_x_F_at_gps_initial_guess(); CHKERRQ(ierr);

	// solving for A_dual with beta's specified
	// ierr = this->solve_grad_part_A_dual(); CHKERRQ(ierr);


	// for (int i=0;i<grid->ndof;i++){

	// 	std::cout << "z values at node i are " << z_grad_A[2*i] << " and " << z_grad_A[2*i+1] << std::endl;
	// }

	// beta's are specified
	// for (int i=0;i<grid->ndof;i++){
	// 	beta[i][0] = 0.0;
	// 	beta[i][1] = 0.0;
	// }

	// beta's are specified
	// ierr = this->solve_curl_part_A_dual(); CHKERRQ(ierr);

	// ierr = this->solve_dual_fields_from_primal(); CHKERRQ(ierr);

	ierr = this->copying_dual_to_petsc_vector(); CHKERRQ(ierr);

	ierr = this->compute_L2_norm_initial_guess(); CHKERRQ(ierr);

	return 0;
}


PetscErrorCode Dual_Solve::copying_dual_to_petsc_vector()
{

	PetscErrorCode ierr;

	if (problem_type == forward){
		throw " copying_dual_to_petsc_vector coded only for 2+4=6 dof (for the inverse problem)";
	}

	int indc[grid->dof_per_node];
	double val[grid->dof_per_node];

	for (int i=0;i<grid->ndof;i++){

		indc[0] = i*grid->dof_per_node;
		indc[1] = i*grid->dof_per_node + 1;

		beta[i][0] = xi_dual[6*i];
		beta[i][1] = xi_dual[6*i+1];

		A_dual[i][0] = xi_dual[6*i+2];
		A_dual[i][1] = xi_dual[6*i+3];
		A_dual[i][2] = xi_dual[6*i+4];
		A_dual[i][3] = xi_dual[6*i+5];

		val[0] = beta[i][0];
		val[1] = beta[i][1];

		for (int k1=0;k1<dim;k1++){
			for (int k2=0;k2<2;k2++){

				val[2+2*k1+k2] = A_dual[i][2*k1+k2];

				indc[2+2*k1+k2] = i*grid->dof_per_node + 2 + 2*k1 + k2;
			}
		}


		ierr = VecSetValues(dual_vec,grid->dof_per_node,indc,val,INSERT_VALUES); CHKERRQ(ierr);

	}// loop over nodes

	// final assembly of global vector (due to multiple proccess)
	ierr = VecAssemblyBegin(dual_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(dual_vec); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode Dual_Solve::copying_petsc_dual_vector_to_stdvector()
{

	PetscErrorCode ierr;

	double max1 = 0.0, max2 = 0.0, max3 = 0.0, max4 = 0.0, max5 = 0.0, max6 = 0.0, max7 = 0.0, max8 = 0.0;

	if (this->count == 0){
		ierr = VecScatterCreateToAll(dual_vec,&ctx_dual,&dual_vec_SEQ); CHKERRQ(ierr);
	}

	ierr = VecScatterBegin(ctx_dual,dual_vec,dual_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_dual,dual_vec,dual_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(dual_vec_SEQ,grid->tdof,ind_set_xi,xi_dual); CHKERRQ(ierr);

	if (problem_type == inverse)
	{
		for (int i=0;i<grid->ndof;i++){

			beta[i][0] = xi_dual[grid->dof_per_node*i];
			beta[i][1] = xi_dual[grid->dof_per_node*i+1];

			A_dual[i][0] = xi_dual[grid->dof_per_node*i+2];
			A_dual[i][1] = xi_dual[grid->dof_per_node*i+3];
			A_dual[i][2] = xi_dual[grid->dof_per_node*i+4];
			A_dual[i][3] = xi_dual[grid->dof_per_node*i+5];

		}// loop over nodes
	}
	else if (problem_type == forward)
	{
		for (int i=0;i<grid->ndof;i++){
			
			beta[i][0] = xi_dual[grid->dof_per_node*i];
			beta[i][1] = xi_dual[grid->dof_per_node*i+1];

			A_dual[i][0] = xi_dual[grid->dof_per_node*i+2];
			A_dual[i][1] = xi_dual[grid->dof_per_node*i+3];
			A_dual[i][2] = xi_dual[grid->dof_per_node*i+4];
			A_dual[i][3] = xi_dual[grid->dof_per_node*i+5];
			A_dual[i][4] = xi_dual[grid->dof_per_node*i+6];
			A_dual[i][5] = xi_dual[grid->dof_per_node*i+7];

		}// loop over nodes
	}

	if (problem_type == inverse)
	{
		for (int i=0;i<grid->ndof;i++){

			xi_dual[grid->dof_per_node*i] = 0.0;
			xi_dual[grid->dof_per_node*i+1] = 0.0;
			xi_dual[grid->dof_per_node*i+2] = 0.0;
			xi_dual[grid->dof_per_node*i+3] = 0.0;
			xi_dual[grid->dof_per_node*i+4] = 0.0;
			xi_dual[grid->dof_per_node*i+5] = 0.0;

		}
	}
	if (problem_type == forward)
	{
		for (int i=0;i<grid->ndof;i++){

			xi_dual[grid->dof_per_node*i] = 0.0;
			xi_dual[grid->dof_per_node*i+1] = 0.0;
			xi_dual[grid->dof_per_node*i+2] = 0.0;
			xi_dual[grid->dof_per_node*i+3] = 0.0;
			xi_dual[grid->dof_per_node*i+4] = 0.0;
			xi_dual[grid->dof_per_node*i+5] = 0.0;
			xi_dual[grid->dof_per_node*i+6] = 0.0;
			xi_dual[grid->dof_per_node*i+7] = 0.0;

		}
	}

	if (this->count == 0){
		ierr = VecScatterCreateToAll(rhs_vec,&ctx_res,&rhs_vec_SEQ); CHKERRQ(ierr);
	}

	ierr = VecScatterBegin(ctx_res,rhs_vec,rhs_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_res,rhs_vec,rhs_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(rhs_vec_SEQ,grid->tdof,ind_set_xi,xi_dual); CHKERRQ(ierr);

	for (int i=0;i<grid->ndof;i++){

		beta_rhs[i][0] = xi_dual[grid->dof_per_node*i];
		beta_rhs[i][1] = xi_dual[grid->dof_per_node*i+1];

		A_dual_rhs[i][0] = xi_dual[grid->dof_per_node*i+2];
		A_dual_rhs[i][1] = xi_dual[grid->dof_per_node*i+3];
		A_dual_rhs[i][2] = xi_dual[grid->dof_per_node*i+4];
		A_dual_rhs[i][3] = xi_dual[grid->dof_per_node*i+5];

		if (problem_type == forward){
			A_dual_rhs[i][4] = xi_dual[grid->dof_per_node*i+6];
			A_dual_rhs[i][5] = xi_dual[grid->dof_per_node*i+7];
		}

		if (std::abs(xi_dual[grid->dof_per_node*i]) > max1){
			max1 = std::abs(xi_dual[grid->dof_per_node*i]);
		}

		if (std::abs(xi_dual[grid->dof_per_node*i + 1]) > max2){
			max2 = std::abs(xi_dual[grid->dof_per_node*i + 1]);
		}

		if (std::abs(xi_dual[grid->dof_per_node*i + 2]) > max3){
			max3 = std::abs(xi_dual[grid->dof_per_node*i + 2]);
		}

		if (std::abs(xi_dual[grid->dof_per_node*i + 3]) > max4){
			max4 = std::abs(xi_dual[grid->dof_per_node*i + 3]);
		}

		if (std::abs(xi_dual[grid->dof_per_node*i + 4]) > max5){
			max5 = std::abs(xi_dual[grid->dof_per_node*i + 4]);
		}

		if (std::abs(xi_dual[grid->dof_per_node*i + 5]) > max6){
			max6 = std::abs(xi_dual[grid->dof_per_node*i + 5]);
		}

		if (problem_type == forward)
		{
			if (std::abs(xi_dual[grid->dof_per_node*i + 6]) > max7){
				max7 = std::abs(xi_dual[grid->dof_per_node*i + 6]);
			}

			if (std::abs(xi_dual[grid->dof_per_node*i + 7]) > max8){
				max8 = std::abs(xi_dual[grid->dof_per_node*i + 7]);
			}
		}

	}// loop over nodes

	std::cout << "norm_rhs_values for " << " beta_1 is " << max1 << " \n "
										<< " beta_2 is " << max2 << " \n "
										<< " Adual_1 is " << max3 << " \n "
										<< " Adual_2 is " << max4 << " \n "
										<< " Adual_3 is " << max5 << " \n "
										<< " Adual_4 is " << max6 << " \n " 
										<< " Adual_5 is " << max7 << " \n " 
										<< " Adual_6 is " << max8 << std::endl;
									


	return 0;
}

PetscErrorCode Dual_Solve::initial_guess_cylindrical()
{


	if (problem_type == forward){
		throw " intial_guess_cylindrical coded only (for the inverse problem)";
	}

	double r = grid->l1;
	double z = grid->l2;

	int ind;

	double r_cr, z_cr; // corrected based on jacobian
	double ratio_cr = 1.0 - lambda_2;
	double area_cr =(PI/2.0)*lambda_1 - 1.0;

	r_cr = r*(1.0 + area_cr); // corrected radii
	z_cr = z*(1.0 - ratio_cr); // corrrected length
	// r_cr = r;
	// z_cr = z;

	std::cout << "r_cr is " << r_cr << std::endl;
	std::cout << "z_cr is " << z_cr << std::endl;

	for(int i=0;i<(grid->nx+1);i++){
		for(int j=0;j<(grid->ny+1);j++){

			// initial data on a plane (manifold)
			ind = i+j*(grid->nx+1);

			x_def[ind][0] = -r_cr + (2.0*r_cr)*((double)i/(double)grid->nx);
			x_def[ind][1] = -z_cr/2.0 + (z_cr)*((double)j/(double)grid->ny);
			// x_def[ind][2] = 0.0;

			// x_def[ind][0] = -PI*r + PI*r*((double)i/(double)grid->nx);
			// x_def[ind][1] = -z/2.0 + (z/2.0)*((double)j/(double)grid->ny);
			// x_def[ind][2] = 0.0;
		}	
	}

	return 0;
}

PetscErrorCode Dual_Solve::initial_guess_plane()
{


	// if (problem_type == forward){
	// 	throw " intial_guess_plane coded only (for the inverse problem)";
	// }

	double r = grid->l1;
	double z = grid->l2;

	int ind;

	double r_cr, z_cr; // corrected based on jacobian
	double area_cr = this->lambda_1*this->lambda_2*(r*z);

	r_cr = r*(this->lambda_1); // corrected radii
	z_cr = area_cr/r_cr; // corrrected length

	std::cout << "r_cr is " << r_cr << std::endl;
	std::cout << "z_cr is " << z_cr << std::endl;

	for(int i=0;i<(grid->nx+1);i++){
		for(int j=0;j<(grid->ny+1);j++){

			// initial data on a plane (manifold)
			ind = i+j*(grid->nx+1);

			x_def[ind][0] = -r_cr/2.0 + r_cr*((double)i/(double)grid->nx);
			x_def[ind][1] = -z_cr/2.0 + z_cr*((double)j/(double)grid->ny);
			
			if (problem_type == forward){
				x_def[ind][2] = 0.0;
			}
			
		}	
	}

	return 0;
}


PetscErrorCode Dual_Solve::initial_guess_trelis()
{

	std::ifstream myfile;

	myfile.open("./input/initial_guess.in",std::ios::in);

	if (problem_type == inverse)
	{
		for(int i=0;i<grid->ndof;i++){

			myfile >> x_def[i][0] >> x_def[i][1];	
		}
	}
	else if (problem_type == forward)
	{
		for(int i=0;i<grid->ndof;i++){

			myfile >> x_def[i][0] >> x_def[i][1] >> x_def[i][2];	
		}
		// scaling the height proportionally
		for (int i=0;i<grid->ndof;++i){
			x_def[i][2] = 0.05*x_def[i][2];
		}
	}

	myfile.close();

	// double r;
	// double r_0; //-0.95
	// double R; 
	// double z_c; 
	// double val;

	// for (int i=0;i<grid->ndof;++i){

	// 	r = 4.5;
	// 	r_0 = -0.5;
	// 	R = std::sqrt(r*r + r_0*r_0);
	// 	z_c = R + ((r_0-R)*2.0*R*(R-r_0)/((R-r_0)*(R-r_0) + r*r));

	// 	val = (R-r_0)*(R-r_0) + x_def[i][0]*x_def[i][0] + x_def[i][1]*x_def[i][1];
	// 	x_def[i][0] = 2.0*R*(R-r_0)*x_def[i][0]/val;//X_ref[i][0];
	// 	x_def[i][1] = 2.0*R*(R-r_0)*x_def[i][1]/val;//X_ref[i][1];
	// 	x_def[i][2] = R + ((r_0-R)*2.0*R*(R-r_0)/val) - z_c;

	// 	// std::cout << " x_def is " << x_def[i][0] << " "
	// 	// 							<< x_def[i][1] << " "
	// 	// 							<< x_def[i][2] << std::endl;
	// }

	return 0;
}


PetscErrorCode Dual_Solve::initial_guess_hat_shape()
{

	double r_cr,r_x,theta,r_def,tmp;
	double r_max = 0.0;

	// lambda1 and lambda2
	double ratio_pr = this->lambda_1/this->lambda_2;

	double area_cr = this->lambda_1*this->lambda_2*undef_global;

	r_def = std::sqrt(area_cr/PI); // for circle as initial guess

	// r_def = std::sqrt(area_cr/(ratio_pr*PI)); // for ellipse as initial guess

	std::cout << "deformed radii is " << r_def << std::endl;

	std::ifstream myfile1;

	myfile1.open("./input/coordinates.in",std::ios::in);

	for (int i=0;i<grid->ndof;++i){

		myfile1 >> x_def[i][0] >> x_def[i][1] >> r_max;
	}

	myfile1.close();

	for (int i=0;i<grid->ndof;++i){

		tmp = std::sqrt(x_def[i][0]*x_def[i][0] +  x_def[i][1]*x_def[i][1]);

		r_max = std::max(r_max,tmp);
	}

	r_cr = r_def/r_max;
	// r_cr = this->lambda_1; // equal lambda1 and lambda2

	std::cout << "r max is " << r_max << " and ratio is " << r_cr << std::endl;

	// scaling to match the jacobian of deformed and undeformed shapes
	for (int i=0;i<grid->ndof;++i){

		r_x = std::sqrt(x_def[i][0]*x_def[i][0] +  x_def[i][1]*x_def[i][1]);

		if (r_x > 1e-10)
		{
			if (x_def[i][1] > 0.0)
				theta = std::acos(x_def[i][0]/r_x);
			else
				theta = 2.0*PI - std::acos(x_def[i][0]/r_x);
		}	
		else
			theta = 0.0;


		x_def[i][0] = r_x*r_cr*std::cos(theta);
		x_def[i][1] = r_x*r_cr*std::sin(theta);

	}

	// for (int i=0;i<grid->ndof;++i){

	// 	x_def[i][0] = (ratio_pr)*x_def[i][0]; // ratio_pr*x
	// 	x_def[i][1] = x_def[i][1]; // y
	// }

	return 0;
}


PetscErrorCode Dual_Solve::get_x_F_at_gps_initial_guess()
{	

	PetscErrorCode ierr;

	MyMat <double> xdef_el(fem->node_per_el,MyVec<double>(this->dim));

	MyMat <double> xdef_0_el(fem->node_per_el,MyVec<double>(this->dim));

	MyMat <double> E_dual_el(fem->ngp,MyVec <double> (2*grid->X_DOF));

	MyMat <double> F_0_dual(dim,MyVec<double>(2));

	MyMat <double> F_0(dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	int ind;

	if (run_type == restart)
	{
		std::ifstream myfile;

		myfile.open("restart_x_0_def.txt",std::ios::in);

		if (problem_type == inverse)
		{
			for(int i=0;i<grid->ndof;i++){
				myfile >> x_0_def[i][0] >> x_0_def[i][1];	
			}
		}
		else if (problem_type == forward)
		{
			for(int i=0;i<grid->ndof;i++){
				myfile >> x_0_def[i][0] >> x_0_def[i][1] >> x_0_def[i][2];	
			}
		}

		myfile.close();
	}

	for (int ie=0;ie<grid->nel;ie++)
	{

		for (int k=0;k<fem->node_per_el;k++){
			for (int d=0;d<this->dim;d++){

				xdef_el[k][d] = x_def[grid->el_conn[ie][k]][d];
				xdef_0_el[k][d] = x_0_def[grid->el_conn[ie][k]][d];

			}
		}

		for (int i=0;i<fem->ngp;i++){
			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
			}
		}


		for (int gp=0;gp<fem->ngp;gp++)
		{

			for (int d=0;d<this->dim;d++){

				x_gp[ie][gp][d] = 0.0; // initializing to zero
				x_0_gp[ie][gp][d] = 0.0; // initializing to zero
				
				for (int k=0;k<fem->ngp;k++){
				
					x_gp[ie][gp][d] = x_gp[ie][gp][d] + fem->psi[gp][k]*xdef_el[k][d];
					x_0_gp[ie][gp][d] = x_0_gp[ie][gp][d] + fem->psi[gp][k]*xdef_0_el[k][d];
				}

				for (int alpha=0;alpha<2;alpha++){
					
					ind = d*2 + alpha;

					F_gp[ie][gp][ind] = 0.0; // initializing to zero
					F_0_gp[ie][gp][ind] = 0.0; // initializing to zero

					for (int k=0;k<fem->ngp;k++){

						F_gp[ie][gp][ind] = F_gp[ie][gp][ind] + fem->dpsi[gp][k][alpha]*xdef_el[k][d]; 
						F_0_gp[ie][gp][ind] = F_0_gp[ie][gp][ind] + fem->dpsi[gp][k][alpha]*xdef_0_el[k][d]; 
					}

					F_0[d][alpha] = F_0_gp[ie][gp][ind];

					// std::cout << "F_0 is " << F_0_gp[ie][gp][ind] << std::endl;

				}

			}

			ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

			ierr = this->tranform_F_to_dual(F_0,G_dual,F_0_dual); CHKERRQ(ierr);

			for (int d=0;d<this->dim;d++){
				for (int alpha=0;alpha<2;alpha++){

					ind = d*2 + alpha;

					F_0_dual_gp[ie][gp][ind] = F_0_dual[d][alpha];
				}
			}


		}// loop over gpts

	}// loop over elems

	return 0;
}


PetscErrorCode Dual_Solve::get_x_F_at_gps_restart_guess()
{	

	PetscErrorCode ierr;

	std::ifstream myfile1, myfile2, myfile3, myfile4;

	myfile1.open("restart_x_gp.txt",std::ios::in);

	if (problem_type == inverse)
	{
		for(int ie=0;ie<grid->nel;ie++){

			for (int gp=0;gp<fem->ngp;gp++){

				myfile1 >> x_gp[ie][gp][0] >> x_gp[ie][gp][1];
			}	
		}
	}
	else if (problem_type == forward)
	{
		for(int ie=0;ie<grid->nel;ie++){

			for (int gp=0;gp<fem->ngp;gp++){

				myfile1 >> x_gp[ie][gp][0] >> x_gp[ie][gp][1] >> x_gp[ie][gp][2];
			}	
		}

	}

	myfile1.close();

	myfile2.open("restart_F_gp.txt",std::ios::in);

	if (problem_type == inverse)
	{
		for(int ie=0;ie<grid->nel;ie++){
			for (int gp=0;gp<fem->ngp;gp++){

				myfile2 >> F_gp[ie][gp][0] >> F_gp[ie][gp][1] >> 
							F_gp[ie][gp][2] >> F_gp[ie][gp][3];
			}	
		}
	}
	else if (problem_type == forward)
	{
		for(int ie=0;ie<grid->nel;ie++){
			for (int gp=0;gp<fem->ngp;gp++){

				myfile2 >> F_gp[ie][gp][0] >> F_gp[ie][gp][1] >> 
							F_gp[ie][gp][2] >> F_gp[ie][gp][3] >>
							F_gp[ie][gp][4] >> F_gp[ie][gp][5];
			}	
		}
	}
	
	myfile2.close();

	myfile3.open("restart_x_0_gp.txt",std::ios::in);

	if (problem_type == inverse)
	{
		for(int ie=0;ie<grid->nel;ie++){
			for (int gp=0;gp<fem->ngp;gp++){

				myfile3 >> x_0_gp[ie][gp][0] >> x_0_gp[ie][gp][1];
			}	
		}
	}
	else if (problem_type == forward)
	{
		for(int ie=0;ie<grid->nel;ie++){
			for (int gp=0;gp<fem->ngp;gp++){

				myfile3 >> x_0_gp[ie][gp][0] >> x_0_gp[ie][gp][1] >> x_0_gp[ie][gp][2];
			}	
		}
	}

	myfile3.close();

	myfile4.open("restart_F_0_gp.txt",std::ios::in);

	if (problem_type == inverse)
	{
		for(int ie=0;ie<grid->nel;ie++){
			for (int gp=0;gp<fem->ngp;gp++){

				myfile4 >> F_0_gp[ie][gp][0] >> F_0_gp[ie][gp][1] >> 
							F_0_gp[ie][gp][2] >> F_0_gp[ie][gp][3];
			}	
		}
	}
	else if (problem_type == forward)
	{
		for(int ie=0;ie<grid->nel;ie++){
			for (int gp=0;gp<fem->ngp;gp++){

				myfile4 >> F_0_gp[ie][gp][0] >> F_0_gp[ie][gp][1] >> 
							F_0_gp[ie][gp][2] >> F_0_gp[ie][gp][3] >>
							F_0_gp[ie][gp][4] >> F_0_gp[ie][gp][5];
			}	
		}
	}
		
	myfile4.close();

	MyMat <double> E_dual_el(fem->ngp,MyVec <double> (2*grid->X_DOF));

	MyMat <double> F_0_dual(this->dim,MyVec<double>(2));

	MyMat <double> F_0(this->dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	int ind;

	for (int ie=0;ie<grid->nel;ie++)
	{

		for (int i=0;i<fem->ngp;i++){
			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
			}
		}


		for (int gp=0;gp<fem->ngp;gp++)
		{

			for (int d=0;d<this->dim;d++){
				
				for (int alpha=0;alpha<2;alpha++){
					
					ind = d*2 + alpha;

					F_0[d][alpha] = F_0_gp[ie][gp][ind];

				}
			}

			ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

			ierr = this->tranform_F_to_dual(F_0,G_dual,F_0_dual); CHKERRQ(ierr);

			for (int d=0;d<this->dim;d++){
				for (int alpha=0;alpha<2;alpha++){

					ind = d*2 + alpha;

					F_0_dual_gp[ie][gp][ind] = F_0_dual[d][alpha];
				}
			}

		}// loop over gpts

	}// loop over elems

	return 0;
}




PetscErrorCode Dual_Solve::deformed_shape_l2_projection()
{

	PetscErrorCode ierr;

	// scalars
	PetscReal rnorm=0.0;
  	PetscInt  its=0;

	int dof_for_xdef_solve = fem->node_per_el*this->dim;

	int indc[dof_for_xdef_solve];

	double **ke;

	ke = new double* [dof_for_xdef_solve];
	for (int i=0;i<dof_for_xdef_solve;i++)
		ke[i]= new double [dof_for_xdef_solve];

	double ke_vec[dof_for_xdef_solve*dof_for_xdef_solve];

	double fe[dof_for_xdef_solve];

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> A_el(fem->node_per_el,MyVec<double>(grid->A_DOF));

	MyMat <double> x_0_gp_el(fem->ngp,MyVec<double>(this->dim));

	ierr = VecSet(F_xdef,0.0);CHKERRQ(ierr);

	int node_el;

	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initialising to zero element stiff matrices and force vectors
		for (int k1=0;k1<dof_for_xdef_solve;k1++){
			for (int k2=0;k2<dof_for_xdef_solve;k2++){

				ke[k1][k2] = 0.0;
			}
			fe[k1] = 0.0;
		}

		// collecting jacobian for the current element
		for (int i=0;i<fem->ngp;i++){			
			detJ_e[i] = fem->detJ_G[ie][i];	
		}

		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<this->dim;k++){
				x_0_gp_el[i][k] = x_0_gp[ie][i][k];
			}
		}

		// A_el for the element
		for (int i=0;i<fem->node_per_el;i++){

			node_el = grid->el_conn[ie][i];

			for (int k=0;k<grid->A_DOF;k++){
				A_el[i][k] = A_dual[node_el][k];
			}
		}


		ierr = this->elem_deformed_geometry_solve_stiff_force(A_el,x_0_gp_el,detJ_e,ke,fe); CHKERRQ(ierr);


		for (int k1=0;k1<fem->node_per_el;k1++){
			for (int d=0;d<this->dim;d++){
				indc[this->dim*k1+d] = this->dim*grid->el_conn[ie][k1] + d;
			}
		}

		// adding values to global matrix and vector
		for (int i=0;i<dof_for_xdef_solve;i++){
			for (int j=0;j<dof_for_xdef_solve;j++){

				ke_vec[i*dof_for_xdef_solve+j] = ke[i][j];
			}
		}

		// adding values to global matrix and vector
		if (this->count == 0){
			ierr = MatSetValues(K_xdef,dof_for_xdef_solve,indc,dof_for_xdef_solve,indc,ke_vec,ADD_VALUES); CHKERRQ(ierr);
		}
		
		ierr = VecSetValues(F_xdef,dof_for_xdef_solve,indc,fe,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_xdef); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_xdef); CHKERRQ(ierr);

	// final assembly of global matrix
	if (this->count == 0){
		ierr = MatAssemblyBegin(K_xdef,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(K_xdef,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}

	// creating the solver context
	ierr = common_utilities->ksp_mumps_solver_petsc(ksp_xdef,K_xdef,pc_xdef,FF_xdef); CHKERRQ(ierr);
	
	//Solve the linear system
	ierr = KSPSolve(ksp_xdef,F_xdef,xdef_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(ksp_xdef,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(ksp_xdef,&its); CHKERRQ(ierr);

	if (::rank==0){
		std::cout << "deformed geometry solve with norm as " << rnorm << " and iterations are " << its << std::endl;
	}

	// context for scattering data to all processors
	if (this->count == 0){
		ierr = VecScatterCreateToAll(xdef_vec,&ctx_xdef,&xdef_vec_SEQ); CHKERRQ(ierr);
	}

	ierr = VecScatterBegin(ctx_xdef,xdef_vec,xdef_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_xdef,xdef_vec,xdef_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(xdef_vec_SEQ,tdof_z,ind_set_z,x_def_vec); CHKERRQ(ierr);



	for (int i=0;i<grid->ndof;i++){
		for (int d=0;d<this->dim;d++){

			x_def[i][d] = x_def_vec[this->dim*i+d];
			// std::cout << " x_def for i " << i << " and d " << d << " is " << x_def[i][d] << std::endl;
		}
	}


	// delete dynamically alotted matrix
	for (int i = 0; i<dof_for_xdef_solve; i++)
		delete[] ke[i];

	delete[] ke;



	return 0;

}

PetscErrorCode Dual_Solve::elem_deformed_geometry_solve_stiff_force(MyMat<double> A_el, MyMat<double> x_0_gp_el,
																	MyVec<double> detJ_e, double **ke, double *fe)
{	

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	int dof_for_xdef_solve = fem->node_per_el*this->dim;

	double K[dof_for_xdef_solve][dof_for_xdef_solve];
	double F[dof_for_xdef_solve];

	int index;

	double tmp;

	// getting x_gp from dual field at all gps
	for (int gp=0;gp<fem->ngp;gp++){	

		// solve for x^k
		for (int k1=0;k1<this->dim;k1++){

			x_gp_el[gp][k1] = 0.0; // initializing to zero
			tmp = 0.0; // initializing to zero

			for (int k2=0;k2<2;k2++){

				index = k1*2 + k2;

				for (int p=0;p<fem->node_per_el;p++){
					tmp = tmp + fem->dpsi[gp][p][k2]*A_el[p][index];
				}
			}
			x_gp_el[gp][k1] = x_0_gp_el[gp][k1] - (1.0/this->px)*tmp;
		} 

	}

	for (int gp=0;gp<fem->ngp;gp++)
	{

		for (int p=0;p<dof_for_xdef_solve;p++){
			for (int q=0;q<dof_for_xdef_solve;q++){
				K[p][q] = 0.0;
			}
			F[p] = 0.0;
		}


		for (int p=0;p<fem->node_per_el;p++){
			for (int q=0;q<fem->node_per_el;q++){

				for (int d=0;d<this->dim;d++)
				{

					K[this->dim*p+d][this->dim*q+d] = fem->psi[gp][p]*fem->psi[gp][q];	

					//K[this->dim*p][this->dim*q+1] = 0.0;
					//K[this->dim*p+1][this->dim*q] = 0.0;
				}// stiffness matrix
				
			}

			for (int d=0;d<this->dim;d++){
				F[this->dim*p + d] = x_gp_el[gp][d]*fem->psi[gp][p];
			}// force vector
			
		}


		for(int p=0;p<dof_for_xdef_solve;p++){
			for(int q=0;q<dof_for_xdef_solve;q++){

				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[gp]*fem->wt[gp];

			}

			fe[p] = fe[p] + F[p]*detJ_e[gp]*fem->wt[gp];
		}


	}// loop over the gauss pts

	return 0;
}

/*
PetscErrorCode Dual_Solve::solve_grad_part_A_dual()
{

	PetscErrorCode ierr;

	// matrices and vector for solving grad part of A
	Mat K_grad;
	Vec F_grad;
	Vec z_vec;

	// Ksp solver, preconditioners
	KSP ksp_z;
	PC pc_z;
	Mat FF_z;

	// scattering data context to sequential vector
	VecScatter 	ctx_z;
	Vec  z_vec_SEQ;

	// scalars
	PetscReal rnorm=0.0;
  	PetscInt  its=0;

	// temporary variables
	int grad_dof = grid->ndof*this->dim;

	ierr = this->mat_create_petsc(K_grad,grad_dof,20,16); CHKERRQ(ierr);

	ierr = this->vec_create_petsc(F_grad,grad_dof); CHKERRQ(ierr);

	ierr = this->vec_create_petsc(z_vec,grad_dof); CHKERRQ(ierr);

	int dof_for_grad_solve = fem->node_per_el*this->dim;

	int indc[dof_for_grad_solve];

	double ke[8][8];

	// ke = new double* [dof_for_grad_solve];
	// for (int i=0;i<dof_for_grad_solve;i++)
	// 	ke[i]= new double [dof_for_grad_solve];

	double fe[dof_for_grad_solve];

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initialising to zero element stiff matrices and force vectors
		for (int k1=0;k1<dof_for_grad_solve;k1++){
			for (int k2=0;k2<dof_for_grad_solve;k2++){

				ke[k1][k2] = 0.0;
			}
			fe[k1] = 0.0;
		}

		// collecting jacobian for the current element
		for (int i=0;i<fem->ngp;i++){			
			detJ_e[i] = fem->detJ_G[ie][i];	
		}

		// storing x_gp fields at gpts for the current element
		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<dim;k++){
				x_gp_el[i][k] = x_gp[ie][i][k];
			}

		}

		ierr = this->elem_grad_A_stiff_force(x_gp_el,detJ_e,ke,fe); CHKERRQ(ierr);


		for (int k1=0;k1<fem->node_per_el;k1++){
			for (int d=0;d<this->dim;d++){
				indc[this->dim*k1+d] = this->dim*grid->el_conn[ie][k1] + d;
			}
		}

		// adding values to global matrix and vector
		ierr = MatSetValues(K_grad,dof_for_grad_solve,indc,dof_for_grad_solve,indc,*ke,ADD_VALUES); CHKERRQ(ierr);
		ierr = VecSetValues(F_grad,dof_for_grad_solve,indc,fe,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_grad); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_grad); CHKERRQ(ierr);

	// final assembly of global matrix
	ierr = MatAssemblyBegin(K_grad,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K_grad,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// boundary conditions
	ierr = this->grad_A_solve_bcs(K_grad,F_grad,z_vec); CHKERRQ(ierr);

	// creating the solver context
	ierr = this->ksp_mumps_solver_petsc(ksp_z,K_grad,pc_z,FF_z); CHKERRQ(ierr);
	
	//Solve the linear system
	ierr = KSPSolve(ksp_z,F_grad,z_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(ksp_z,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(ksp_z,&its); CHKERRQ(ierr);

	if (::rank==0){
		std::cout << "grad part (z) of A solve with norm as " << rnorm << " and iterations are " << its << std::endl;
	}

	// context for scattering data to all processors
	ierr = VecScatterCreateToAll(z_vec,&ctx_z,&z_vec_SEQ); CHKERRQ(ierr);

	ierr = VecScatterBegin(ctx_z,z_vec,z_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_z,z_vec,z_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(z_vec_SEQ,tdof_z,ind_set_z,z_grad_A); CHKERRQ(ierr);


	// delete objects
	MatDestroy(&K_grad);	VecDestroy(&F_grad); VecDestroy(&z_vec);

	KSPDestroy(&ksp_z); VecScatterDestroy(&ctx_z); VecDestroy(&z_vec_SEQ);

	// delete dynamically alotted matrix
	// for (int i = 0; i<dof_for_grad_solve; i++)
	// 	delete[] ke[i];

	// delete[] ke;



	return 0;
}
*/

/*
PetscErrorCode Dual_Solve::elem_grad_A_stiff_force(MyMat<double> x_gp_el, MyVec<double> detJ_e,
													double ke[][8], double *fe)
{	

	int dof_for_grad_solve = fem->node_per_el*this->dim;

	double K[dof_for_grad_solve][dof_for_grad_solve];
	double F[dof_for_grad_solve];

	for (int gp=0;gp<fem->ngp;gp++)
	{

		for (int p=0;p<dof_for_grad_solve;p++){
			for (int q=0;q<dof_for_grad_solve;q++){
				K[p][q] = 0.0;
			}
			F[p] = 0.0;
		}


		for (int p=0;p<fem->node_per_el;p++){
			for (int q=0;q<fem->node_per_el;q++){

				for (int d=0;d<this->dim;d++)
				{

					K[this->dim*p+d][this->dim*q+d] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0] 
															+ fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];	

					//K[this->dim*p][this->dim*q+1] = 0.0;
					//K[this->dim*p+1][this->dim*q] = 0.0;
				}// stiffness matrix
				
			}

			for (int d=0;d<this->dim;d++){
				F[this->dim*p + d] = this->px*x_gp_el[gp][d];
			}// force vector
			
		}


		for(int p=0;p<dof_for_grad_solve;p++){
			for(int q=0;q<dof_for_grad_solve;q++){

				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[gp]*fem->wt[gp];

			}

			fe[p] = fe[p] + F[p]*detJ_e[gp]*fem->wt[gp];
		}


	}// loop over the gauss pts

	return 0;
}
*/

/*
PetscErrorCode Dual_Solve::grad_A_solve_bcs(Mat &K_grad, Vec &F_grad, Vec &z_vec)
{

	PetscErrorCode ierr;

	int tbel = 3;
  	int rows[tbel];
  	double val[tbel];

  	rows[0] = 2*0;
  	rows[1] = 2*0+1;
  	rows[2] = 2*grid->nx+1;

  	val[0] = 0.0;
  	val[1] = 0.0;
  	val[2] = 0.0;

	ierr = VecSetValues(z_vec,tbel,rows,val,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValues(F_grad,tbel,rows,val,INSERT_VALUES); CHKERRQ(ierr);

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_grad); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_grad); CHKERRQ(ierr);

	// final assembly of global vector
	ierr = VecAssemblyBegin(z_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(z_vec); CHKERRQ(ierr);

	
	ierr = MatZeroRows(K_grad,tbel,rows,1.0,z_vec,F_grad); CHKERRQ(ierr);

	return 0;
}
*/

/*
PetscErrorCode Dual_Solve::solve_curl_part_A_dual()
{	

	PetscErrorCode ierr;

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

	// scalars
	PetscReal rnorm=0.0;
  	PetscInt  its=0;

	// temporary variables
	int xi_dof = grid->ndof*grid->A_DOF;

	ierr = this->mat_create_petsc(K_xi,xi_dof,40,30); CHKERRQ(ierr);

	ierr = this->vec_create_petsc(F_xi,xi_dof); CHKERRQ(ierr);

	ierr = this->vec_create_petsc(xi_vec,xi_dof); CHKERRQ(ierr);

	int dof_for_xi_solve = fem->node_per_el*grid->A_DOF;

	int indc[dof_for_xi_solve];

	double ke[16][16];

	int index;

	// ke = new double* [dof_for_xi_solve];
	// for (int i=0;i<dof_for_xi_solve;i++)
	// 	ke[i]= new double [dof_for_xi_solve];

	double fe[dof_for_xi_solve];

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	MyMat <double> F_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat<double> beta_el(fem->node_per_el,MyVec<double> (grid->BETA_DOF));

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	// loop over the elements
	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initialising to zero element stiff matrices and force vectors
		for (int k1=0;k1<dof_for_xi_solve;k1++){
			for (int k2=0;k2<dof_for_xi_solve;k2++){

				ke[k1][k2] = 0.0;
			}
			fe[k1] = 0.0;
		}

		// collecting jacobian for the current element
		for (int i=0;i<fem->ngp;i++){			
			
			detJ_e[i] = fem->detJ_G[ie][i];

			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
				//std::cout << "E_dual_el for i and k is " << E_dual_el[i][k] << std::endl;
			}	
		}

		// beta_el is specified
		for (int i=0;i<fem->node_per_el;i++){
			
			for (int k=0;k<grid->BETA_DOF;k++){
				beta_el[i][k] = beta[grid->el_conn[ie][i]][k];
			}

		}


		// storing x_gp and F_gp fields at gpts for the current element
		for (int i=0;i<fem->ngp;i++){

			for (int k1=0;k1<this->dim;k1++){
				for (int k2=0;k2<2;k2++){

					index = 2*k1 + k2;

					F_gp_el[i][index] = F_gp[ie][i][index];	
				}
			}

			for (int k=0;k<this->dim;k++){
				x_gp_el[i][k] = x_gp[ie][i][k];
			}

		}

		// ierr = this->elem_curl_A_stiff_force(F_gp_el,detJ_e,E_dual_el,beta_el,ke,fe); CHKERRQ(ierr);

		ierr = this->elem_L2_A_stiff_force(x_gp_el,F_gp_el,detJ_e,E_dual_el,beta_el,ke,fe); CHKERRQ(ierr);

		for (int i=0;i<fem->node_per_el;i++){
			for (int k1=0;k1<this->dim;k1++){
				for (int k2=0;k2<2;k2++){

					index = k1*2 + k2;
					indc[F_DOF*i+index] = F_DOF*grid->el_conn[ie][i] + index;	
				}
			}
		}

		// adding values to global matrix and vector
		ierr = MatSetValues(K_xi,dof_for_xi_solve,indc,dof_for_xi_solve,indc,*ke,ADD_VALUES); CHKERRQ(ierr);
		ierr = VecSetValues(F_xi,dof_for_xi_solve,indc,fe,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_xi); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_xi); CHKERRQ(ierr);

	// final assembly of global matrix
	ierr = MatAssemblyBegin(K_xi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K_xi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// boundary conditions
	ierr = this->A_solve_bcs(K_xi,F_xi,xi_vec); CHKERRQ(ierr);

	// creating the solver context
	ierr = this->ksp_mumps_solver_petsc(ksp_xi,K_xi,pc_xi,FF_xi); CHKERRQ(ierr);
	
	//Solve the linear system
	ierr = KSPSolve(ksp_xi,F_xi,xi_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(ksp_xi,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(ksp_xi,&its); CHKERRQ(ierr);

	if (::rank==0){
		std::cout << "A solve with norm as " << rnorm << " and iterations are " << its << std::endl;
	}

	// context for scattering data to all processors
	ierr = VecScatterCreateToAll(xi_vec,&ctx_xi,&xi_vec_SEQ); CHKERRQ(ierr);

	ierr = VecScatterBegin(ctx_xi,xi_vec,xi_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_xi,xi_vec,xi_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(xi_vec_SEQ,tdof_xi,ind_set_xi,xi_curl_A); CHKERRQ(ierr);


	// delete objects
	MatDestroy(&K_xi);	VecDestroy(&F_xi); VecDestroy(&xi_vec);

	KSPDestroy(&ksp_xi); VecScatterDestroy(&ctx_xi); VecDestroy(&xi_vec_SEQ);

	// // delete dynamically alotted matrix
	// for (int i = 0; i<dof_for_xi_solve; i++)
	// 	delete[] ke[i];

	// delete[] ke;


	return 0;

}
*/

/*
PetscErrorCode Dual_Solve::solve_dual_fields_from_primal()
{	

	PetscErrorCode ierr;

	// scalars
	PetscReal rnorm=0.0;
  	PetscInt  its=0;

  	double fdoteigv = 0.0;

	// temporary variables
	int xi_dof = grid->ndof*grid->dof_per_node;

	ierr = this->mat_create_petsc(K_xi,xi_dof,60,50); CHKERRQ(ierr);

	ierr = this->vec_create_petsc(F_xi,xi_dof); CHKERRQ(ierr);

	ierr = this->vec_create_petsc(xi_vec,xi_dof); CHKERRQ(ierr);

	ierr = VecSet(xi_vec,0.0); CHKERRQ(ierr);

	int dof_for_xi_solve = fem->node_per_el*grid->dof_per_node;

	int indc[dof_for_xi_solve];

	double **ke; //[24][24]

	int index;

	ke = new double* [dof_for_xi_solve];
	for (int i=0;i<dof_for_xi_solve;i++)
		ke[i]= new double [dof_for_xi_solve];

	double fe[dof_for_xi_solve];

	double ke_vec[dof_for_xi_solve*dof_for_xi_solve];

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	MyMat <double> F_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> x_0_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> F_0_dual_gp_el(fem->ngp,MyVec<double>(F_DOF));

	// loop over the elements
	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initialising to zero element stiff matrices and force vectors
		for (int k1=0;k1<dof_for_xi_solve;k1++){
			for (int k2=0;k2<dof_for_xi_solve;k2++){

				ke[k1][k2] = 0.0;
			}
			fe[k1] = 0.0;
		}

		// collecting jacobian for the current element
		for (int i=0;i<fem->ngp;i++){			
			
			detJ_e[i] = fem->detJ_G[ie][i];

			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
				//std::cout << "E_dual_el for i and k is " << E_dual_el[i][k] << std::endl;
			}	
		}


		// storing x_gp and F_gp fields at gpts for the current element
		for (int i=0;i<fem->ngp;i++){

			for (int k1=0;k1<this->dim;k1++){
				for (int k2=0;k2<2;k2++){

					index = 2*k1 + k2;

					F_gp_el[i][index] = F_gp[ie][i][index];	
				}
			}

			for (int k1=0;k1<this->dim;k1++){
				for (int k2=0;k2<2;k2++){

					index = k1*2 + k2;
					
					F_0_dual_gp_el[i][index] = F_0_dual_gp[ie][i][index];
				}
			}

			for (int k=0;k<this->dim;k++){
				x_gp_el[i][k] = x_gp[ie][i][k];
			}

			for (int k=0;k<this->dim;k++){
				x_0_gp_el[i][k] = x_0_gp[ie][i][k];
			}

		}

		
		// ierr = this->elem_curl_A_stiff_force(F_gp_el,detJ_e,E_dual_el,beta_el,ke,fe); CHKERRQ(ierr);

		// ierr = this->elem_L2_A_stiff_force(x_gp_el,F_gp_el,detJ_e,E_dual_el,beta_el,ke,fe); CHKERRQ(ierr);

		ierr = this->elem_L2_dual_stiff_force(x_gp_el,x_0_gp_el,F_gp_el,F_0_dual_gp_el,detJ_e,E_dual_el,ke,fe); CHKERRQ(ierr);

		for (int i=0;i<fem->node_per_el;i++){
			for (int k1=0;k1<grid->dof_per_node;k1++){

				indc[grid->dof_per_node*i+k1] = grid->dof_per_node*grid->el_conn[ie][i] + k1;	
				
			}
		}

		// adding values to global matrix and vector
		for (int i=0;i<dof_for_xi_solve;i++){
			for (int j=0;j<dof_for_xi_solve;j++){

				ke_vec[i*dof_for_xi_solve+j] = ke[i][j];
			}
		}

		ierr = MatSetValues(K_xi,dof_for_xi_solve,indc,dof_for_xi_solve,indc,ke_vec,ADD_VALUES); CHKERRQ(ierr);
		ierr = VecSetValues(F_xi,dof_for_xi_solve,indc,fe,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_xi); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_xi); CHKERRQ(ierr);

	// final assembly of global matrix
	ierr = MatAssemblyBegin(K_xi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K_xi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// PetscViewer writer1;

	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"stiffness_matrix.txt",&writer1); CHKERRQ(ierr);
	// ierr = PetscViewerPushFormat(writer1,PETSC_VIEWER_ASCII_INDEX); CHKERRQ(ierr);
	// ierr = MatView(K_xi,writer1); CHKERRQ(ierr);

	// linear solve using petsc
	// {
	// 	// boundary conditions
	// 	ierr = this->dual_solve_bcs(K_xi,F_xi,xi_vec); CHKERRQ(ierr);

	// 	// creating the solver context
	// 	ierr = this->ksp_mumps_solver_petsc(ksp_xi,K_xi,pc_xi,FF_xi); CHKERRQ(ierr);
		
	// 	//Solve the linear system
	// 	ierr = KSPSolve(ksp_xi,F_xi,xi_vec);CHKERRQ(ierr);

	// 	ierr = KSPGetResidualNorm(ksp_xi,&rnorm); CHKERRQ(ierr);
	// 	ierr = KSPGetTotalIterations(ksp_xi,&its); CHKERRQ(ierr);

	// }

	// bcs on x
	// ierr = this->dual_solve_bcs(K_xi,F_xi,xi_vec); CHKERRQ(ierr);

	int nevls = 550; //grid->tdof

	// eigenvalue solver
	{	

		ierr = MatCreateVecs(K_xi,NULL,&eigvAr);
		ierr = MatCreateVecs(K_xi,NULL,&eigvAi);

		ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

		ierr = EPSSetOperators(eps,K_xi,NULL); CHKERRQ(ierr);

		ierr = EPSSetProblemType(eps,EPS_NHEP); CHKERRQ(ierr);

		// EPSGetRG(eps,&rg);
		// RGSetType(rg,RGINTERVAL);
		// RGIntervalSetEndpoints(rg,-100.0,100.0,-0.01,0.01);
		// RGSetComplement(rg,PETSC_TRUE);
		// sort eigenvalue approximations wrt a target, otherwise convergence will be erratic 
		// EPSSetTarget(eps,1.0);
		// EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);

		ierr = EPSSetDimensions(eps,nevls,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
		ierr = EPSSetTolerances(eps,1e-7,2000); CHKERRQ(ierr); 

		ierr = EPSSetFromOptions(eps);

		ierr = EPSSolve(eps);

		if (::rank==0){
			std::cout << "dual solve with for eigenvalues done " << std::endl;
		}

	}

	for (int i=0;i<nevls;i++){

		ierr = EPSGetEigenpair(eps,i,&kr,&ki,eigvAr,eigvAi); CHKERRQ(ierr);

		ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error); CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"eigen value is %9E+ %9Ei , %12g\n",(double)kr,(double)ki,(double)error); CHKERRQ(ierr);

		fdoteigv = 0.0;

		ierr = VecDot(eigvAr,F_xi,&fdoteigv); CHKERRQ(ierr);

		fdoteigv = fdoteigv/(double)kr;

		ierr = VecAXPY(xi_vec,fdoteigv,eigvAr); CHKERRQ(ierr);

	}


	// context for scattering data to all processors
	ierr = VecScatterCreateToAll(xi_vec,&ctx_xi,&xi_vec_SEQ); CHKERRQ(ierr);

	ierr = VecScatterBegin(ctx_xi,xi_vec,xi_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_xi,xi_vec,xi_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(xi_vec_SEQ,tdof_xi,ind_set_xi,xi_dual); CHKERRQ(ierr);

	// deleting objects
	ierr = EPSDestroy(&eps); CHKERRQ(ierr);

	ierr = VecDestroy(&eigvAr); CHKERRQ(ierr);

	ierr = VecDestroy(&eigvAi); CHKERRQ(ierr);

	ierr = this->delete_petsc_objects(); CHKERRQ(ierr);
	
	// delete dynamically alotted matrix
	for (int i = 0; i<dof_for_xi_solve; i++)
		delete[] ke[i];

	delete[] ke;


	return 0;
	
}
*/


PetscErrorCode Dual_Solve::compute_L2_norm_initial_guess()
{


	PetscErrorCode ierr;

	MyMat<double> beta_el(fem->node_per_el,MyVec<double> (grid->BETA_DOF));

	MyMat<double> A_el(fem->node_per_el, MyVec<double>(grid->A_DOF));

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> x_0_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> F_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> F_0_dual_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	double norm_l2_x = 0.0; // mapping equation for x

	double norm_l2_F = 0.0; // mapping equation for F

	int index;

	for (int ie=0;ie<grid->nel;++ie)
	{	

		// getting dual fields for the nodes of each elem 
		for (int i=0;i<fem->node_per_el;i++){
			
			for (int k=0;k<grid->BETA_DOF;k++){
				beta_el[i][k] = beta[grid->el_conn[ie][i]][k];
			}

			for (int k=0;k<grid->A_DOF;k++){
				A_el[i][k] = A_dual[grid->el_conn[ie][i]][k];
			}
		}

		// getting basis vectors and jacobian for each elem
		for (int i=0;i<fem->ngp;i++){
			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
			}
			
			detJ_e[i] = fem->detJ_G[ie][i];	
		}

		// getting x_gp and F_gp fields at gpts for the current element
		for (int i=0;i<fem->ngp;i++){

			for (int k1=0;k1<this->dim;k1++){
				for (int k2=0;k2<2;k2++){

					index = 2*k1 + k2;

					F_gp_el[i][index] = F_gp[ie][i][index];	
				}
			}


			for (int k1=0;k1<this->dim;k1++){
				for (int k2=0;k2<2;k2++){

					index = 2*k1 + k2;

					F_0_dual_gp_el[i][index] = F_0_dual_gp[ie][i][index];	
				}
			}

			for (int k=0;k<this->dim;k++){
				x_gp_el[i][k] = x_gp[ie][i][k];
			}

			for (int k=0;k<this->dim;k++){
				x_0_gp_el[i][k] = x_0_gp[ie][i][k];
			}

		}

		ierr = this->elem_compute_L2_norm_initial_guess(x_gp_el,x_0_gp_el,F_gp_el,F_0_dual_gp_el,beta_el,A_el,E_dual_el,
															detJ_e,&norm_l2_x,&norm_l2_F); CHKERRQ(ierr);


		norm_global_x[ie] = norm_l2_x;

		norm_global_F[ie] = norm_l2_F;


	}// loop over elems 


	return 0;
}


PetscErrorCode Dual_Solve::elem_compute_L2_norm_initial_guess(MyMat<double> x_gp_el, MyMat<double> x_0_gp_el, 
															MyMat<double> F_gp_el, MyMat<double> F_0_dual_gp_el,
															MyMat<double> beta_el,MyMat<double> A_el, MyMat<double> E_dual_el,
															MyVec<double> detJ_e,double *norm_l2_x, double *norm_l2_F)

{
	if (problem_type == forward){
		throw " elem_compute_L2_norm_initial_guess coded only (for the inverse problem)";
	}

	PetscErrorCode ierr;

	*norm_l2_x = 0.0;

	*norm_l2_F = 0.0;

	double norm_x_gp, norm_F_gp;	

	int index;

	MyMat <double> F_dual(dim,MyVec<double>(2));

	MyMat <double> F_real(dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	MyMat <double> G_real(2,MyVec<double>(2));

	MyMat <double> C(2,MyVec<double>(2));

	MyMat <double> C_inv(2,MyVec<double>(2));

	MyMat <double> d_mu1_F(fem->ngp,MyVec<double>(2*this->dim));

	MyMat <double> d_mu2_F(fem->ngp,MyVec<double>(2*this->dim));

	MyMat <double> A_rhs(fem->ngp,MyVec<double>(2*this->dim));

	double A_mat_gp[2*this->dim];

	double div_A[this->dim];

	double A_at_gp;

	double beta_1_gp, beta_2_gp;

	double tr_C = 0.0, det_C = 0.0;

	double dmu1_dF = 0.0;

	double dmu2_dF = 0.0;

	// rhs at all gpts
	for (int gp=0;gp<fem->ngp;gp++)
	{

		for (int k=0;k<this->dim;k++){
			for (int gamma=0;gamma<2;gamma++){

				index = k*2 + gamma;

				F_real[k][gamma] = F_gp_el[gp][index];
			}
		}


		ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);

		ierr = this->tranform_F_to_dual(F_real,G_dual,F_dual); CHKERRQ(ierr);

		ierr = this->get_C(F_real,F_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F_real,G_real,&det_C); CHKERRQ(ierr);

		// ierr = this->tranform_C_on_both_legs(G_real,G_dual,C,C_tranf); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(C,C_inv); CHKERRQ(ierr); // in tranformed basis

		for (int k1=0; k1<this->dim; k1++){
			for (int k2=0; k2<2; k2++){

				index = k1*2 + k2;

				ierr = this->get_derivative_of_mus(C_inv,tr_C,det_C,G_dual,F_real,k1,k2,&dmu1_dF,&dmu2_dF); CHKERRQ(ierr);

				A_rhs[gp][index] = this->pf*(F_dual[k1][k2] - F_0_dual_gp_el[gp][index]);

				d_mu1_F[gp][index] = dmu1_dF;

				d_mu2_F[gp][index] = dmu2_dF;

			}
		}

	} // rhs at all gpts


	// loop over the gpts
	for (int gp=0;gp<fem->ngp;gp++){


		// norms
		norm_x_gp = 0.0;

		norm_F_gp = 0.0;

		// A at gp
		for (int k1=0; k1<dim; k1++){
			for (int k2=0; k2<2; k2++){

				index = k1*2 + k2;

				ierr = this->get_A_field_at_gp(A_el,index,gp,&A_at_gp); CHKERRQ(ierr);

				A_mat_gp[index] = A_at_gp; 
			}
		}


		// beta at gp
		ierr = this->get_beta_field_at_gp(beta_el,gp,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);	

		// div of A at gp
		ierr = this->get_div_A_field_at_gp(A_el,gp,div_A); CHKERRQ(ierr);
							

		for (int i=0;i<this->dim;i++){

			norm_x_gp = norm_x_gp + (div_A[i] + this->px*(x_gp_el[gp][i] - x_0_gp_el[gp][i]))*
											(div_A[i] + this->px*(x_gp_el[gp][i] - x_0_gp_el[gp][i]));
		}	


		for (int i=0;i<2*this->dim;i++){

			norm_F_gp = norm_F_gp + (A_mat_gp[i] + d_mu1_F[gp][i]*beta_1_gp + d_mu2_F[gp][i]*beta_2_gp + A_rhs[gp][i])*
									(A_mat_gp[i] + d_mu1_F[gp][i]*beta_1_gp + d_mu2_F[gp][i]*beta_2_gp + A_rhs[gp][i]);
		}			 

		*norm_l2_x = *norm_l2_x + norm_x_gp; // *detJ_e[gp]*fem->wt[gp]
		*norm_l2_F = *norm_l2_F + norm_F_gp; // *detJ_e[gp]*fem->wt[gp]	


	} // loop over the gpts


	return 0;
}

// PetscErrorCode Dual_Solve::elem_curl_A_stiff_force(MyMat <double> F_gp_el,MyVec<double> detJ_e,
// 													MyMat <double> E_dual_el, MyMat<double> beta_el,
// 													double ke[][16], double *fe)
// {

// 	PetscErrorCode ierr;

// 	int dof_for_xi_solve = fem->node_per_el*grid->A_DOF;

// 	double K[dof_for_xi_solve][dof_for_xi_solve];

// 	double F[dof_for_xi_solve];

// 	MyMat<double> G_dual(2,MyVec<double>(2));	

// 	int index;

// 	MyVec <double> curl_A_tilde(2);

// 	MyMat <double> F_dual(dim,MyVec<double>(2));

// 	MyMat <double> F_real(dim,MyVec<double>(2));

// 	MyMat <double> G_dual(2,MyVec<double>(2));

// 	MyMat <double> G_real(2,MyVec<double>(2));

// 	MyMat <double> C(2,MyVec<double>(2));

// 	MyMat <double> C_tranf(2,MyVec<double>(2));

// 	MyMat <double> C_inv(2,MyVec<double>(2));

// 	double tr_C = 0.0, det_C = 0.0;

// 	double d_mu1_F = 0.0;

// 	double d_mu2_F = 0.0;

// 	double beta_1_gp = 0.0;

// 	double beta_2_gp = 0.0;

// 	MyMat <double> A_rhs(fem->ngp,MyVec<double>(this->dim*2));

// 	// rhs at all gpts
// 	for (int gp=0;gp<fem->ngp;gp++){

// 		for (int k=0;k<this->dim;k++){
// 			for (int gamma=0;gamma<2;gamma++){

// 				index = k*2 + gamma;

// 				F_real[k][gamma] = F_gp_el[gp][index];
// 			}
// 		}


// 		ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

// 		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);

// 		ierr = this->tranform_F_to_dual(F_real,G_dual,F_dual); CHKERRQ(ierr);

// 		ierr = this->get_C(F_real,F_dual,C); CHKERRQ(ierr);

// 		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

// 		ierr = this->determinant_of_C(F_real,G_real,&det_C); CHKERRQ(ierr);

// 		ierr = this->tranform_C_on_both_legs(G_real,G_dual,C,C_tranf); CHKERRQ(ierr);

// 		ierr = fem->my_mat_inverse(C,C_inv); CHKERRQ(ierr); // in tranformed basis

// 		ierr = this->get_beta_field_at_gp(beta_el,gp,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);

// 		for (int k1=0; k1<this->dim; k1++){
// 			for (int k2=0; k2<2; k2++){

// 				index = k1*2 + k2;

// 				ierr = this->get_derivative_of_mus(C_inv,tr_C,det_C,G_dual,F_real,k1,k2,&d_mu1_F,&d_mu2_F); CHKERRQ(ierr);

// 				A_rhs[gp][index] = this->pf*F_dual[k1][k2] + d_mu1_F*beta_1_gp + d_mu2_F*beta_2_gp;

// 			}
// 		}

// 	} // rhs at all gpts


// 	// curl A_tilde at the center of domain
// 	for (int i=0;i<fem->ngp;i++){

// 		curl_A_tilde[0] = curl_A_tilde[0] + (A_rhs[i][1]*fem->dpsi_c[i][0])
// 										 	- (A_rhs[i][0]*fem->dpsi_c[i][1]);

// 		curl_A_tilde[1] = curl_A_tilde[1] + (A_rhs[i][3]*fem->dpsi_c[i][0])
// 										 	- (A_rhs[i][2]*fem->dpsi_c[i][1]);
// 	}

// 	for (int gp=0;gp<fem->ngp;gp++){

// 		for (int k1=0;k1<dof_for_xi_solve;k1++){
// 			for (int k2=0;k2<dof_for_xi_solve;k2++){

// 				K[k1][k2] = 0.0;
// 			}

// 			F[k1] = 0.0;
// 		}


// 		// elem stiff matrix and force vector
// 		for (int p=0;p<grid->node_per_el;p++){
// 			for (int q=0;q<grid->node_per_el;q++){

// 				K[F_DOF*p][F_DOF*q] = fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0];
// 				K[F_DOF*p][F_DOF*q+1] = -fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1];
// 				K[F_DOF*p][F_DOF*q+2] = 0.0;
// 				K[F_DOF*p][F_DOF*q+3] = 0.0;

// 				K[F_DOF*p+1][F_DOF*q] = -fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0];
// 				K[F_DOF*p+1][F_DOF*q+1] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];
// 				K[F_DOF*p+1][F_DOF*q+2] = 0.0;
// 				K[F_DOF*p+1][F_DOF*q+3] = 0.0;

// 				K[F_DOF*p+2][F_DOF*q] = 0.0;
// 				K[F_DOF*p+2][F_DOF*q+1] = 0.0;
// 				K[F_DOF*p+2][F_DOF*q+2] = fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0];
// 				K[F_DOF*p+2][F_DOF*q+3] = -fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1];

// 				K[F_DOF*p+3][F_DOF*q] = 0.0;
// 				K[F_DOF*p+3][F_DOF*q+1] = 0.0;
// 				K[F_DOF*p+3][F_DOF*q+2] = -fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1] + fem->dpsi[gp][p][1]*dem->dpsi[gp][q][0];
// 				K[F_DOF*p+3][F_DOF*q+3] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];
// 			}

// 			F[F_DOF*p] = curl_A_tilde[0]*fem->dpsi[gp][p][1];
// 			F[F_DOF*p+1] = -curl_A_tilde[0]*fem->dpsi[gp][p][0];
// 			F[F_DOF*p+2] = curl_A_tilde[1]*fem->dpsi[gp][p][1];
// 			F[F_DOF*p+3] = -curl_A_tilde[1]*fem->dpsi[gp][p][0];
// 		}


// 		for (int p=0;p<dof_for_xi_solve;p++){
// 			for (int q=0;q<dof_for_xi_solve;q++){
				
// 				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[gp]*fem->wt[gp];		
// 			}

// 			fe[p] = fe[p] + F[p]*detJ_e[gp]*fem->wt[gp];
// 		}


// 	} // loop over the gpts


// 	return 0;

// }


/*
PetscErrorCode Dual_Solve::elem_L2_A_stiff_force(MyMat <double> x_gp_el,
												MyMat <double> F_gp_el, MyVec<double> detJ_e,
												MyMat <double> E_dual_el, MyMat<double> beta_el,
												double ke[][16], double *fe)
{

	PetscErrorCode ierr;

	int dof_for_xi_solve = fem->node_per_el*grid->A_DOF;

	double K[dof_for_xi_solve][dof_for_xi_solve];

	double F[dof_for_xi_solve];	

	int index;

	MyMat <double> F_dual(dim,MyVec<double>(2));

	MyMat <double> F_real(dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	MyMat <double> G_real(2,MyVec<double>(2));

	MyMat <double> C(2,MyVec<double>(2));

	MyMat <double> C_tranf(2,MyVec<double>(2));

	MyMat <double> C_inv(2,MyVec<double>(2));

	double tr_C = 0.0, det_C = 0.0;

	double d_mu1_F = 0.0;

	double d_mu2_F = 0.0;

	double beta_1_gp = 0.0;

	double beta_2_gp = 0.0;

	MyMat <double> A_rhs(fem->ngp,MyVec<double>(this->dim*2));

	// rhs at all gpts
	for (int gp=0;gp<fem->ngp;gp++){

		for (int k=0;k<this->dim;k++){
			for (int gamma=0;gamma<2;gamma++){

				index = k*2 + gamma;

				F_real[k][gamma] = F_gp_el[gp][index];
			}
		}


		ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);

		ierr = this->tranform_F_to_dual(F_real,G_dual,F_dual); CHKERRQ(ierr);

		ierr = this->get_C(F_real,F_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F_real,G_real,&det_C); CHKERRQ(ierr);

		// ierr = this->tranform_C_on_both_legs(G_real,G_dual,C,C_tranf); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(C,C_inv); CHKERRQ(ierr); // in tranformed basis

		ierr = this->get_beta_field_at_gp(beta_el,gp,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);

		for (int k1=0; k1<this->dim; k1++){
			for (int k2=0; k2<2; k2++){

				index = k1*2 + k2;

				ierr = this->get_derivative_of_mus(C_inv,tr_C,det_C,G_dual,F_real,k1,k2,&d_mu1_F,&d_mu2_F); CHKERRQ(ierr);

				A_rhs[gp][index] = this->pf*F_dual[k1][k2] + d_mu1_F*beta_1_gp + d_mu2_F*beta_2_gp;

				// std::cout << "F_dual and F_real for k1 and k2 " << k1 
				// 					<< " and " << k2 << " is " << F_dual[k1][k2] << " and " << F_real[k1][k2] << std::endl;

				// std::cout << "A_rhs for k1 and k2 " << k1 
				// 					<< " and " << k2 << " is " << A_rhs[gp][index] << std::endl;
			}
		}

	} // rhs at all gpts


	// loop over the gpts
	for (int gp=0;gp<fem->ngp;gp++){

		for (int k1=0;k1<dof_for_xi_solve;k1++){
			for (int k2=0;k2<dof_for_xi_solve;k2++){

				K[k1][k2] = 0.0;
			}

			F[k1] = 0.0;
		}


		// elem stiff matrix and force vector
		for (int p=0;p<fem->node_per_el;p++){
			for (int q=0;q<fem->node_per_el;q++){

				K[F_DOF*p][F_DOF*q] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0];
				K[F_DOF*p][F_DOF*q+1] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1];
				K[F_DOF*p][F_DOF*q+2] = 0.0;
				K[F_DOF*p][F_DOF*q+3] = 0.0;

				K[F_DOF*p+1][F_DOF*q] = fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0];
				K[F_DOF*p+1][F_DOF*q+1] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];
				K[F_DOF*p+1][F_DOF*q+2] = 0.0;
				K[F_DOF*p+1][F_DOF*q+3] = 0.0;

				K[F_DOF*p+2][F_DOF*q] = 0.0;
				K[F_DOF*p+2][F_DOF*q+1] = 0.0;
				K[F_DOF*p+2][F_DOF*q+2] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0];
				K[F_DOF*p+2][F_DOF*q+3] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1];

				K[F_DOF*p+3][F_DOF*q] = 0.0;
				K[F_DOF*p+3][F_DOF*q+1] = 0.0;
				K[F_DOF*p+3][F_DOF*q+2] = fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0];
				K[F_DOF*p+3][F_DOF*q+3] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];

				// std::cout << "K 0 and 0 entry is " << K[0][0] << std::endl;
			}

			F[F_DOF*p] = -(A_rhs[gp][0]*fem->psi[gp][p] + this->px*x_gp_el[gp][0]*fem->dpsi[gp][p][0]);
			F[F_DOF*p+1] = -(A_rhs[gp][1]*fem->psi[gp][p] + this->px*x_gp_el[gp][0]*fem->dpsi[gp][p][1]);
			F[F_DOF*p+2] = -(A_rhs[gp][2]*fem->psi[gp][p] + this->px*x_gp_el[gp][1]*fem->dpsi[gp][p][0]);
			F[F_DOF*p+3] = -(A_rhs[gp][3]*fem->psi[gp][p] + this->px*x_gp_el[gp][1]*fem->dpsi[gp][p][1]);

			// std::cout << "F entry is " << F[F_DOF*p]  << " , "
			// 							 << F[F_DOF*p+1]  << " , "
			// 							 << F[F_DOF*p+2]  << " , "
			// 							 << F[F_DOF*p+3]  << std::endl;
		}


		for (int p=0;p<dof_for_xi_solve;p++){
			for (int q=0;q<dof_for_xi_solve;q++){
				
				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[gp]*fem->wt[gp];		
			}

			fe[p] = fe[p] + F[p]*detJ_e[gp]*fem->wt[gp];
		}


	} // loop ove the gpts


	return 0;
}*/

/*
PetscErrorCode Dual_Solve::A_solve_bcs(Mat &K_xi, Vec &F_xi, Vec &xi_vec)
{

	PetscErrorCode ierr;

	int tbel = 4*grid->nbel;
  	int rows[4];
  	double val[4];

  	int rows_full[tbel];

  	for (int i=0;i<grid->nbel;i++){

  		rows[0] = 4*grid->b_el_conn[i][0];
	  	rows[1] = 4*grid->b_el_conn[i][0]+1;
	  	rows[2] = 4*grid->b_el_conn[i][0]+2;
	  	rows[3] = 4*grid->b_el_conn[i][0]+3;

	  	val[0] = 0.0;
	  	val[1] = 0.0;
	  	val[2] = 0.0;
	  	val[3] = 0.0;

	  	rows_full[4*i] = rows[0];
	  	rows_full[4*i+1] = rows[1];
	  	rows_full[4*i+2] = rows[2];
	  	rows_full[4*i+3] = rows[3];

		ierr = VecSetValues(xi_vec,4,rows,val,INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValues(F_xi,4,rows,val,INSERT_VALUES); CHKERRQ(ierr);
	}

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_xi); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_xi); CHKERRQ(ierr);

	// final assembly of global vector
	ierr = VecAssemblyBegin(xi_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xi_vec); CHKERRQ(ierr);
	
	ierr = MatZeroRows(K_xi,tbel,rows_full,1.0,xi_vec,F_xi); CHKERRQ(ierr);

	return 0;
}
*/

PetscErrorCode Dual_Solve::elem_L2_dual_stiff_force(MyMat <double> x_gp_el, MyMat <double> x_0_gp_el,
												MyMat <double> F_gp_el, MyMat <double> F_0_dual_gp_el, MyVec<double> detJ_e,
												MyMat <double> E_dual_el, double **ke, double *fe)
{

	if (problem_type == forward){
		throw " least square stiffness matrix and force vector coded only (for the inverse problem)";
	}

	PetscErrorCode ierr;

	int dof_for_xi_solve = fem->node_per_el*grid->dof_per_node;

	double K[dof_for_xi_solve][dof_for_xi_solve];

	double F[dof_for_xi_solve];	

	int index;

	int T_DOF = grid->dof_per_node;

	MyMat <double> F_dual(dim,MyVec<double>(2));

	MyMat <double> F_real(dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	MyMat <double> G_real(2,MyVec<double>(2));

	MyMat <double> C(2,MyVec<double>(2));

	// MyMat <double> C_tranf(2,MyVec<double>(2));

	MyMat <double> C_inv(2,MyVec<double>(2));

	MyMat <double> d_mu1_F(fem->ngp,MyVec<double>(2*this->dim));

	MyMat <double> d_mu2_F(fem->ngp,MyVec<double>(2*this->dim));

	double tr_C = 0.0, det_C = 0.0;

	double dmu1_dF = 0.0;

	double dmu2_dF = 0.0;

	double coeff_b1_b1, coeff_b1_b2, coeff_b2_b2;

	MyMat <double> A_rhs(fem->ngp,MyVec<double>(this->dim*2));

	// rhs at all gpts
	for (int gp=0;gp<fem->ngp;gp++)
	{

		for (int k=0;k<this->dim;k++){
			for (int gamma=0;gamma<2;gamma++){

				index = k*2 + gamma;

				F_real[k][gamma] = F_gp_el[gp][index];
			}
		}


		ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

		// std::cout << "G_dual for gp " << gp << " is " << G_dual[0][0] << ", "
		// 						<< G_dual[0][1] << ", "
		// 						<< G_dual[1][0] << ", "
		// 						<< G_dual[1][1] << std::endl;

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);

		ierr = this->tranform_F_to_dual(F_real,G_dual,F_dual); CHKERRQ(ierr);

		ierr = this->get_C(F_real,F_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F_real,G_real,&det_C); CHKERRQ(ierr);

		// ierr = this->tranform_C_on_both_legs(G_real,G_dual,C,C_tranf); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(C,C_inv); CHKERRQ(ierr); // in tranformed basis

		for (int k1=0; k1<this->dim; k1++){
			for (int k2=0; k2<2; k2++){

				index = k1*2 + k2;

				ierr = this->get_derivative_of_mus(C_inv,tr_C,det_C,G_dual,F_real,k1,k2,&dmu1_dF,&dmu2_dF); CHKERRQ(ierr);

				A_rhs[gp][index] = this->pf*(F_dual[k1][k2] - F_0_dual_gp_el[gp][index]);

				d_mu1_F[gp][index] = dmu1_dF;

				d_mu2_F[gp][index] = dmu2_dF;

			}
		}

	} // rhs at all gpts


	// loop over the gpts
	for (int gp=0;gp<fem->ngp;gp++){

		for (int k1=0;k1<dof_for_xi_solve;k1++){
			for (int k2=0;k2<dof_for_xi_solve;k2++){

				K[k1][k2] = 0.0;
			}

			F[k1] = 0.0;
		}

		coeff_b1_b1 =  d_mu1_F[gp][0]*d_mu1_F[gp][0] + d_mu1_F[gp][1]*d_mu1_F[gp][1] +
							d_mu1_F[gp][2]*d_mu1_F[gp][2] + d_mu1_F[gp][3]*d_mu1_F[gp][3];

		coeff_b2_b2 =  d_mu2_F[gp][0]*d_mu2_F[gp][0] + d_mu2_F[gp][1]*d_mu2_F[gp][1] +
							d_mu2_F[gp][2]*d_mu2_F[gp][2] + d_mu2_F[gp][3]*d_mu2_F[gp][3];

		coeff_b1_b2 =  d_mu1_F[gp][0]*d_mu2_F[gp][0] + d_mu1_F[gp][1]*d_mu2_F[gp][1] +
							d_mu1_F[gp][2]*d_mu2_F[gp][2] + d_mu1_F[gp][3]*d_mu2_F[gp][3];


		// elem stiff matrix and force vector
		for (int p=0;p<fem->node_per_el;p++){
			for (int q=0;q<fem->node_per_el;q++){

				
				K[T_DOF*p][T_DOF*q] = coeff_b1_b1*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p][T_DOF*q+1] = coeff_b1_b2*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p][T_DOF*q+2] = d_mu1_F[gp][0]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p][T_DOF*q+3] = d_mu1_F[gp][1]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p][T_DOF*q+4] = d_mu1_F[gp][2]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p][T_DOF*q+5] = d_mu1_F[gp][3]*fem->psi[gp][p]*fem->psi[gp][q];


				K[T_DOF*p+1][T_DOF*q] = coeff_b1_b2*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+1][T_DOF*q+1] = coeff_b2_b2*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+1][T_DOF*q+2] = d_mu2_F[gp][0]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+1][T_DOF*q+3] = d_mu2_F[gp][1]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+1][T_DOF*q+4] = d_mu2_F[gp][2]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+1][T_DOF*q+5] = d_mu2_F[gp][3]*fem->psi[gp][p]*fem->psi[gp][q];


				K[T_DOF*p+2][T_DOF*q] = d_mu1_F[gp][0]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+2][T_DOF*q+1] = d_mu2_F[gp][0]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+2][T_DOF*q+2] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0];
				K[T_DOF*p+2][T_DOF*q+3] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1];
				K[T_DOF*p+2][T_DOF*q+4] = 0.0;
				K[T_DOF*p+2][T_DOF*q+5] = 0.0;

				K[T_DOF*p+3][T_DOF*q] = d_mu1_F[gp][1]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+3][T_DOF*q+1] = d_mu2_F[gp][1]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+3][T_DOF*q+2] = fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0];
				K[T_DOF*p+3][T_DOF*q+3] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];
				K[T_DOF*p+3][T_DOF*q+4] = 0.0;
				K[T_DOF*p+3][T_DOF*q+5] = 0.0;

				K[T_DOF*p+4][T_DOF*q] = d_mu1_F[gp][2]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+4][T_DOF*q+1] = d_mu2_F[gp][2]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+4][T_DOF*q+2] = 0.0;
				K[T_DOF*p+4][T_DOF*q+3] = 0.0;
				K[T_DOF*p+4][T_DOF*q+4] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][0]*fem->dpsi[gp][q][0];
				K[T_DOF*p+4][T_DOF*q+5] = fem->dpsi[gp][p][0]*fem->dpsi[gp][q][1];

				K[T_DOF*p+5][T_DOF*q] = d_mu1_F[gp][3]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+5][T_DOF*q+1] = d_mu2_F[gp][3]*fem->psi[gp][p]*fem->psi[gp][q];
				K[T_DOF*p+5][T_DOF*q+2] = 0.0;
				K[T_DOF*p+5][T_DOF*q+3] = 0.0;
				K[T_DOF*p+5][T_DOF*q+4] = fem->dpsi[gp][p][1]*fem->dpsi[gp][q][0];
				K[T_DOF*p+5][T_DOF*q+5] = fem->psi[gp][p]*fem->psi[gp][q] + fem->dpsi[gp][p][1]*fem->dpsi[gp][q][1];


			}

			F[T_DOF*p] = -(A_rhs[gp][0]*d_mu1_F[gp][0] + A_rhs[gp][1]*d_mu1_F[gp][1] +
							A_rhs[gp][2]*d_mu1_F[gp][2] + A_rhs[gp][3]*d_mu1_F[gp][3])*fem->psi[gp][p];
			
			F[T_DOF*p+1] = -(A_rhs[gp][0]*d_mu2_F[gp][0] + A_rhs[gp][1]*d_mu2_F[gp][1] +
							A_rhs[gp][2]*d_mu2_F[gp][2] + A_rhs[gp][3]*d_mu2_F[gp][3])*fem->psi[gp][p];

			F[T_DOF*p+2] = -(A_rhs[gp][0]*fem->psi[gp][p] + this->px*(x_gp_el[gp][0] - x_0_gp_el[gp][0])*fem->dpsi[gp][p][0]);
			F[T_DOF*p+3] = -(A_rhs[gp][1]*fem->psi[gp][p] + this->px*(x_gp_el[gp][0] - x_0_gp_el[gp][0])*fem->dpsi[gp][p][1]);
			F[T_DOF*p+4] = -(A_rhs[gp][2]*fem->psi[gp][p] + this->px*(x_gp_el[gp][1] - x_0_gp_el[gp][1])*fem->dpsi[gp][p][0]);
			F[T_DOF*p+5] = -(A_rhs[gp][3]*fem->psi[gp][p] + this->px*(x_gp_el[gp][1] - x_0_gp_el[gp][1])*fem->dpsi[gp][p][1]);


		}


		for (int p=0;p<dof_for_xi_solve;p++){
			for (int q=0;q<dof_for_xi_solve;q++){
				
				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[gp]*fem->wt[gp];		
			}

			fe[p] = fe[p] + F[p]*detJ_e[gp]*fem->wt[gp];
		}


	} // loop over the gpts


	return 0;
}


PetscErrorCode Dual_Solve::dual_solve_bcs(Mat &K_xi, Vec &F_xi, Vec &xi_vec)
{

	PetscErrorCode ierr;

	int tbel = grid->nbel*2;
  	int rows[2];
  	double val[2];

  	int rows_full[tbel];

  	for (int i=0;i<grid->nbel;i++){

  		rows[0] = grid->dof_per_node*grid->b_el_conn[i][0];
  		rows[1] = grid->dof_per_node*grid->b_el_conn[i][0]+1;
  		// rows[2] = grid->dof_per_node*grid->b_el_conn[i][0]+2;
	  	// rows[3] = grid->dof_per_node*grid->b_el_conn[i][0]+3;
	  	// rows[4] = grid->dof_per_node*grid->b_el_conn[i][0]+4;
	  	// rows[5] = grid->dof_per_node*grid->b_el_conn[i][0]+5;

	  	val[0] = 0.0;
	  	val[1] = 0.0;
	  	// val[2] = 0.0;
	  	// val[3] = 0.0;
	  	// val[4] = 0.0;
	  	// val[5] = 0.0;

	  	rows_full[2*i] = rows[0];
	  	rows_full[2*i+1] = rows[1];
	  	// rows_full[6*i+2] = rows[2];
	  	// rows_full[6*i+3] = rows[3];
	  	// rows_full[6*i+4] = rows[4];
	  	// rows_full[6*i+5] = rows[5];

		ierr = VecSetValues(xi_vec,2,rows,val,INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValues(F_xi,2,rows,val,INSERT_VALUES); CHKERRQ(ierr);
	}

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_xi); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_xi); CHKERRQ(ierr);

	// final assembly of global vector
	ierr = VecAssemblyBegin(xi_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xi_vec); CHKERRQ(ierr);
	
	ierr = MatZeroRows(K_xi,tbel,rows_full,1.0,xi_vec,F_xi); CHKERRQ(ierr);

	return 0;
}




PetscErrorCode Dual_Solve::boundary_conditions()
{	

	PetscErrorCode ierr;

	double vals[4];
	int rows[4];

	for (int ibe=0;ibe<grid->nbel;ibe++){

		rows[0] = grid->dof_b_el_conn[ibe][2];
		rows[1] = grid->dof_b_el_conn[ibe][3];
		rows[2] = grid->dof_b_el_conn[ibe][4];
		rows[3] = grid->dof_b_el_conn[ibe][5];

		vals[0] = 0.0;
		vals[1] = 0.0;
		vals[2] = 0.0;
		vals[3] = 0.0;

		ierr = VecSetValues(rhs_vec,4,rows,vals,INSERT_VALUES); CHKERRQ(ierr);

	}

	ierr = VecAssemblyBegin(rhs_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs_vec); CHKERRQ(ierr);

	return 0;
}


PetscErrorCode Dual_Solve::force_vector_assembly()
{

	PetscErrorCode ierr;

	ierr = VecSet(rhs_vec,0.0); CHKERRQ(ierr);

	if (this->count == 0){
		ierr = VecSet(mass_vec,0.0); CHKERRQ(ierr);
	}

	MyMat<double> beta_el(fem->node_per_el,MyVec<double> (grid->BETA_DOF));

	MyMat<double> A_el(fem->node_per_el, MyVec<double>(grid->A_DOF));

	//MyVec <int> node_el(fem->node_per_el);

	int node_el;

	double rhs_el[fem->dof_per_el];

	double mass_el[fem->dof_per_el];

	int indc[fem->dof_per_el];

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> x_0_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> F_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> F_guess_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> F_0_dual_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	MyMat <double> xdef_el(fem->node_per_el,MyVec<double>(this->dim));

	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initializing rhs_el and mass_el to be zero (passed as pointers: should not sum)
		for (int i=0;i<fem->dof_per_el;i++){
			rhs_el[i] = 0.0;
			mass_el[i] = 0.0;
		}

		// getting dual fields for each elem to get mapping 
		for (int i=0;i<fem->node_per_el;i++){

			node_el = grid->el_conn[ie][i];
			
			for (int k=0;k<grid->BETA_DOF;k++){
				beta_el[i][k] = beta[node_el][k];
			}

			for (int k=0;k<grid->A_DOF;k++){
				A_el[i][k] = A_dual[node_el][k];
			}

			for (int k=0;k<this->dim;k++){
				xdef_el[i][k] = x_def[node_el][k];
			}
		}

		// getting basis vectors and jacobian for each elem
		for (int i=0;i<fem->ngp;i++){
			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
				//std::cout << "E_dual_el for i and k is " << E_dual_el[i][k] << std::endl;
			}
			
			detJ_e[i] = fem->detJ_G[ie][i];	
		}

		// std::cout<< "********************" << std::endl;
		// std::cout<< "element no is " << ie << std::endl;
		// std::cout<< "********************" << std::endl;

		// getting indicies for assembly of elem vector to global vector
		for (int k=0;k<(grid->A_DOF + grid->BETA_DOF)*fem->node_per_el;k++){
			indc[k] = grid->dof_el_conn[ie][k];
		}

		// storing primal fields at gpts from previous iteration as guess and base state values (fixed)
		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<F_DOF;k++){
				F_guess_gp_el[i][k] = F_gp[ie][i][k];
			}

			for (int k=0;k<F_DOF;k++){
				F_0_dual_gp_el[i][k] = F_0_dual_gp[ie][i][k];
			}

			for (int k=0;k<this->dim;k++){
				x_0_gp_el[i][k] = x_0_gp[ie][i][k];
			}
		}


		// primal fields (obtained from mapping from dual to primal)
		ierr = this->get_primal_fields(beta_el,A_el,E_dual_el,xdef_el,F_guess_gp_el,
												x_0_gp_el,x_gp_el,F_0_dual_gp_el,F_gp_el); CHKERRQ(ierr);

		// rhs vector per elem
		ierr = this->elem_force_vector(beta_el,A_el,E_dual_el,detJ_e,rhs_el,x_gp_el,F_gp_el); CHKERRQ(ierr);


		if (this->count == 0){
			// mass vector per element
			ierr = this->elem_mass_vector(detJ_e,mass_el); CHKERRQ(ierr);
		}

		// storing primal fields at gpts
		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<dim;k++){
				x_gp[ie][i][k] = x_gp_el[i][k];
			}

			for (int k=0;k<F_DOF;k++){
				F_gp[ie][i][k] = F_gp_el[i][k];
			}
		}	

		// assembly to global rhs vector
		ierr = VecSetValues(rhs_vec,fem->dof_per_el,indc,rhs_el,ADD_VALUES); CHKERRQ(ierr);

		if (this->count == 0){
			// assembly to global mass vector
			ierr = VecSetValues(mass_vec,fem->dof_per_el,indc,mass_el,ADD_VALUES); CHKERRQ(ierr);
		}


	} // loop over the elements

	// final assembly of global vector (due to multiple proccess)
	ierr = VecAssemblyBegin(rhs_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs_vec); CHKERRQ(ierr);

	if (this->count == 0){
		// final assembly of global vector (due to multiple proccess)
		ierr = VecAssemblyBegin(mass_vec); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(mass_vec); CHKERRQ(ierr);
	}

	return 0;

}


PetscErrorCode Dual_Solve::elem_mass_vector(MyVec<double> detJ_e, double* mass_el)
{

	double Me[fem->dof_per_el];

	for (int i=0;i<fem->ngp;i++)
	{		

		for (int k=0;k<fem->dof_per_el;k++){
			Me[k] = 0.0;
		}

		for (int p=0;p<fem->node_per_el;++p){

			for (int q=0;q<grid->dof_per_node;++q){

				Me[p*grid->dof_per_node+q] = fem->psi[i][p];

			}
		}

		for(int p=0;p<fem->dof_per_el;p++){

			mass_el[p] = mass_el[p] + Me[p]*detJ_e[i]*fem->wt[i];
		}

	}// loop over gpts

	return 0;

}

PetscErrorCode Dual_Solve::elem_force_vector(MyMat <double> beta_el, MyMat <double> A_el, 
											 MyMat <double> E_dual_el,  MyVec<double> detJ_e, double* rhs_el, 
											 MyMat <double> x_gp_el,  MyMat <double> F_gp_el)
{

	PetscErrorCode ierr;

	MyMat<double> F(dim,MyVec<double>(2));

	MyMat<double> F_dual(dim,MyVec<double>(2));

	MyMat<double> G_dual(2,MyVec<double>(2));

	MyMat<double> G_real(2,MyVec<double>(2));

	MyMat<double> C(2,MyVec<double>(2));

	double tr_C, det_C;

	double mu1, mu2;

	int tmp;

	double fe[fem->dof_per_el];

	for (int i=0;i<fem->ngp;i++)
	{	

		for (int k1=0;k1<dim;k1++){
			for(int k2=0;k2<2;k2++){

				tmp = k1*2 + k2;
				F[k1][k2] = F_gp_el[i][tmp];
			}
		}

		for (int k=0;k<fem->dof_per_el;k++){
			fe[k] = 0.0;
		}

		ierr = this->get_G_dual(i,E_dual_el,G_dual); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);		

		ierr = this->tranform_F_to_dual(F,G_dual,F_dual); CHKERRQ(ierr);

		ierr = this->get_C(F,F_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F,G_real,&det_C); CHKERRQ(ierr);
		
		ierr = this->get_invariants(tr_C,det_C,&mu1,&mu2); CHKERRQ(ierr);

		// std::cout << "mu1 is " << mu1 << " and mu2 is " << mu2 << std::endl;
		// std::cout << "mu1 imposed is " << lambda_1*lambda_1 << " and mu2 imposed is " << lambda_2*lambda_2 << std::endl; 

		for (int p=0;p<fem->node_per_el;p++){

			fe[p*grid->dof_per_node] = fem->psi[i][p]*(mu1- (lambda_1*lambda_1));
			fe[p*grid->dof_per_node+1] = fem->psi[i][p]*(mu2- (lambda_2*lambda_2));

			for (int k1=0;k1<dim;k1++){
				for (int k2=0;k2<2;k2++){

					fe[p*grid->dof_per_node+2+k1*2+k2] = fem->psi[i][p]*F[k1][k2] 
														+ fem->dpsi[i][p][k2]*x_gp_el[i][k1]; 
				}
			}
		}

		for (int p=0;p<fem->dof_per_el;p++){

			rhs_el[p] = rhs_el[p] + fe[p]*detJ_e[i]*fem->wt[i];
		}

	}// loop over gpts

	// std::cout << "rhs_el 0 dof is " <<  rhs_el[0] << std::endl;
	// std::cout << "rhs_el 1 dof is " <<  rhs_el[1] << std::endl;
	// std::cout << "rhs_el 2 dof is " <<  rhs_el[2] << std::endl;
	// std::cout << "rhs_el 3 dof is " <<  rhs_el[3] << std::endl;
	// std::cout << "rhs_el 4 dof is " <<  rhs_el[4] << std::endl;
	// std::cout << "rhs_el 5 dof is " <<  rhs_el[5] << std::endl;
	// std::cout << "rhs_el 6 dof is " <<  rhs_el[6] << std::endl;
	// std::cout << "rhs_el 7 dof is " <<  rhs_el[7] << std::endl;

	return 0;
}

PetscErrorCode Dual_Solve::get_primal_fields(MyMat <double> beta_el, MyMat <double> A_el,
											MyMat <double> E_dual_el, MyMat <double> xdef_el,
											MyMat <double> F_guess_gp_el, 
											MyMat <double> x_0_gp_el, MyMat <double> & x_gp_el,
											MyMat <double> F_0_dual_gp_el, MyMat <double> & F_gp_el)
{

	PetscErrorCode ierr;

	int index;

	int solution;

	MyMat<double> F_sol(dim,MyVec<double>(2));

	MyMat<double> F_guess(dim,MyVec<double>(2));

	MyMat<double> F_0_dual(dim,MyVec<double>(2)); // base state in dual

	double tmp;

	// for (int i=0;i<fem->ngp;i++){
	// 	for (int k=0;k<2*grid->X_DOF;k++){
	// 		std::cout << "E_dual_el for i= " << i << " and k= " << k << " is " << E_dual_el[i][k] << std::endl;
	// 	}
	// }

	for (int gp=0;gp<fem->ngp;gp++){	

		// solve for x^k
		for (int k1=0;k1<this->dim;k1++){

			x_gp_el[gp][k1] = 0.0; // initializing to zero
			tmp = 0.0;

			for (int k2=0;k2<2;k2++){

				index = k1*2 + k2;

				for (int p=0;p<fem->node_per_el;p++){
					tmp = tmp + fem->dpsi[gp][p][k2]*A_el[p][index]; 
				}
			}
			x_gp_el[gp][k1] = x_0_gp_el[gp][k1] - (1.0/this->px)*tmp;
			
			// std::cout << "diff in x and x_0 " << - (1.0/this->px)*tmp << std::endl;
		} 		

		// for (int d=0;d<this->dim;d++){

		// 	x_gp_el[gp][d] = 0.0; // initializing to zero
			
		// 	for (int k=0;k<fem->ngp;k++){
			
		// 		x_gp_el[gp][d] = x_gp_el[gp][d] + fem->psi[gp][k]*xdef_el[k][d];
		// 	}

		// }

		// for (int k1=0;k1<this->dim;k1++){
		// 	for (int k2=0; k2<2; k2++){

		// 		F_guess[k1][k2] = 0.0;

		// 		for (int p=0;p<fem->node_per_el;p++){
		// 			F_guess[k1][k2] = F_guess[k1][k2] + fem->dpsi[gp][p][k2]*xdef_el[p][k1];
		// 		}
		// 	}
		// }


		for (int k1=0;k1<this->dim;k1++){
			for (int k2=0;k2<2;k2++){
				
				index = k1*2 + k2;
				F_guess[k1][k2] = F_guess_gp_el[gp][index];
				F_0_dual[k1][k2] = F_0_dual_gp_el[gp][index];
			}
		}

		// solve for F
		// ierr = this->picards_iteration_for_F(gp,beta_el,A_el,E_dual_el,F_guess,F_0_dual,F_sol); CHKERRQ(ierr);

		ierr = this->local_newton_for_F(gp,beta_el,A_el,E_dual_el,F_guess,F_0_dual,F_sol,&solution); CHKERRQ(ierr);
		// std::cout << "value of solution is " << solution << std::endl;

		// if (solution == 0){

		// 	for (int k1=0;k1<this->dim;k1++){
		// 		for (int k2=0;k2<2;k2++){
					
		// 			index = k1*2 + k2; 
		// 			if (k1 == 0 && k2 == 0)
		// 				F_guess[k1][k2] = F_guess_gp_el[gp][index] + 0.5;
		// 			else if (k1 == 0 && k2 == 1)
		// 				F_guess[k1][k2] = F_guess_gp_el[gp][index] + 0.1;
		// 			else if (k1 == 1 && k2 == 0)
		// 				F_guess[k1][k2] = F_guess_gp_el[gp][index] + 1.5;
		// 			else if (k1 == 1 && k2 == 1)
		// 				F_guess[k1][k2] = F_guess_gp_el[gp][index] + 0.4;

		// 		}
		// 	}

		// 	ierr = this->local_newton_for_F(gp,beta_el,A_el,E_dual_el,F_guess,F_0_dual,F_sol,&solution); CHKERRQ(ierr);
		// }

		// if (solution == 0){
		// 	ierr = this->picards_iteration_for_F(gp,beta_el,A_el,E_dual_el,F_guess,F_0_dual,F_sol); CHKERRQ(ierr);
		// }

		// updating the solution
		for (int k1=0;k1<dim;k1++){
			for (int k2=0;k2<2;k2++){

				index = k1*2 + k2;
				F_gp_el[gp][index] = F_sol[k1][k2];	
			}
		}

		// for (int k1=0;k1<dim;k1++){
		// 	for (int k2=0;k2<2;k2++){

		// 		index = k1*2 + k2;
		// 		std::cout << " F diff after mapping " << F_sol[k1][k2] - F_guess[k1][k2] << std::endl;
		// 	}
		// }	

	} // loop over gpts



	return 0;
}


PetscErrorCode Dual_Solve::picards_iteration_for_F(int gp, MyMat <double> beta_el, MyMat <double> A_el,
													MyMat <double> E_dual_el, MyMat <double> F_guess,
													MyMat <double> F_0_dual, MyMat <double> &F_sol)
{

	PetscErrorCode ierr;

	double residual_local = 1.0;

	double delta_F_norm = 1.0;

	int index;

	int local_step = 0;

	double A_at_gp = 0.0;

	double d_mu1_F = 0.0, d_mu2_F = 0.0;

	double beta_1_gp = 0.0, beta_2_gp = 0.0;

	MyMat <double> rhs_local(this->dim,MyVec<double>(2));

	MyMat <double> delta_F(this->dim,MyVec<double>(2));

	MyMat <double> F_guess_dual(dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	MyMat <double> G_real(2,MyVec<double>(2));

	MyMat <double> C(2,MyVec<double>(2));

	MyMat <double> C_tranf(2,MyVec<double>(2));

	MyMat <double> C_inv(2,MyVec<double>(2));

	double tr_C = 0.0, det_C = 0.0;

	// double h_term = 0.0;

	double mu1 = 0.0, mu2 = 0.0;

	ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

	ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);

	ierr = this->get_beta_field_at_gp(beta_el,gp,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);

	ierr = this->tranform_F_to_dual(F_guess,G_dual,F_guess_dual); CHKERRQ(ierr);

	while (delta_F_norm > tol_local || residual_local > tol_local)
	{	

		if (local_step > 0){
			ierr = this->tranform_F_from_dual(F_guess_dual,G_real,F_guess); CHKERRQ(ierr);
		}

		ierr = this->get_C(F_guess,F_guess_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F_guess,G_real,&det_C); CHKERRQ(ierr);

		ierr = this->get_invariants(tr_C,det_C,&mu1,&mu2); CHKERRQ(ierr);

		// ierr = this->tranform_C_on_both_legs(G_real,G_dual,C,C_tranf); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(C,C_inv); CHKERRQ(ierr); // in tranformed basis




		for (int k1=0; k1<dim; k1++){
			for (int k2=0; k2<2; k2++){

				index = k1*2 + k2;

				ierr = this->get_derivative_of_mus(C_inv,tr_C,det_C,G_dual,F_guess,k1,k2,&d_mu1_F,&d_mu2_F); CHKERRQ(ierr);

				// ierr = this->get_H_term(F_guess,G_dual,k1,k2,&h_term); CHKERRQ(ierr);

				ierr = this->get_A_field_at_gp(A_el,index,gp,&A_at_gp); CHKERRQ(ierr);

				rhs_local[k1][k2] = (1.0/this->pf)*(A_at_gp  
											+ d_mu1_F*beta_1_gp
											+ d_mu2_F*beta_2_gp
											+ this->pf*(F_guess_dual[k1][k2] - F_0_dual[k1][k2])); 

				// if (residual_local > 1e+3){
				// 	std::cout << "derivatives of mu for k1 and k2 " << k1 << " and " << k2 << " is " << d_mu1_F << " and " << d_mu2_F << std::endl;
				// 	std::cout << "A for k1 and k2 " << k1 << " and " << k2 << " is " << A_at_gp << std::endl;
				// 	std::cout << "beta_1_gp is " << beta_1_gp << std::endl;
				// 	std::cout << "beta_2_gp is " << beta_2_gp << std::endl;
				// }
				
				delta_F[k1][k2] = -(1.0/this->pf)*(A_at_gp + d_mu1_F*beta_1_gp + d_mu2_F*beta_2_gp)
											 - (F_guess_dual[k1][k2] - F_0_dual[k1][k2]); 										 
			}
		}

		// delta_F_norm and residual_local
		delta_F_norm = 0.0; 
		residual_local = 0.0;

		for (int k1=0; k1<dim; k1++){
			for (int k2=0; k2<2; k2++){
				 delta_F_norm = delta_F_norm +  delta_F[k1][k2]*delta_F[k1][k2];
				 residual_local =  residual_local + rhs_local[k1][k2]*rhs_local[k1][k2];
			}
		}

		delta_F_norm = std::sqrt(delta_F_norm);
		residual_local = std::sqrt(residual_local);

		// new guess
		for (int k1=0; k1<dim; k1++){
			for (int k2=0; k2<2; k2++){

				F_guess_dual[k1][k2] = F_guess_dual[k1][k2] + delta_F[k1][k2];
				 
			}
		}

		local_step = local_step + 1;

		if (residual_local > 1e+3 || local_step > 500){
		std::cout << "********************" << std::endl;
		std::cout << " local step for gp = " << gp << " is " << local_step << std::endl;
		std::cout << " local norm is " << delta_F_norm << " and " << residual_local << std::endl;
		std::cout << "********************" << std::endl;

		}

		if (local_step > 1000){
			std::cout << "local norms are " << delta_F_norm  << " and " << residual_local << std::endl;

			throw "local_step for gp convergence is greater than 1000";
		}


	} // while loop ends

	// tranform back to required basis (to use in gradient flow)
	ierr = this->tranform_F_from_dual(F_guess_dual,G_real,F_sol); CHKERRQ(ierr);

	return 0;

}

PetscErrorCode Dual_Solve::local_newton_for_F(int gp, MyMat <double> beta_el, MyMat <double> A_el,
												MyMat <double> E_dual_el, MyMat <double> F_guess,
												MyMat <double> F_0_dual, MyMat <double> &F_sol, int *solution)
{

	PetscErrorCode ierr;

	//**** here the F is solved in real basis unline the 
	//**** Picard's iteration. Hence, no transformation is
	//**** required when solution is obtained. 

	*solution = 1;

	MyMat <double> k_local(F_DOF,MyVec <double>(F_DOF));

	MyVec <double> f_local(F_DOF);

	double residual_local = 1.0;

	double delta_F_norm = 1.0;

	int index1, index2;

	int local_step = 0;

	double A_at_gp = 0.0;

	MyVec <double> d_mu1_F(F_DOF);

	MyVec <double> d_mu2_F(F_DOF);

	MyMat <double> dd_mu1_FF(F_DOF,MyVec <double> (F_DOF));

	MyMat <double> dd_mu2_FF(F_DOF,MyVec <double> (F_DOF));

	double beta_1_gp = 0.0, beta_2_gp = 0.0;

	MyVec <double> delta_F(F_DOF);

	MyMat <double> F_guess_dual(dim,MyVec<double>(2));

	MyMat <double> F_0_real(dim,MyVec<double>(2));

	MyMat <double> G_dual(2,MyVec<double>(2));

	MyMat <double> G_real(2,MyVec<double>(2));

	MyMat <double> C(2,MyVec<double>(2));

	MyMat <double> C_tranf(2,MyVec<double>(2));

	MyMat <double> C_inv(2,MyVec<double>(2));

	MyVec <double> eig_v1(2);

	MyVec <double> eig_v2(2);

	double tr_C = 0.0, det_C = 0.0;

	// double h_term = 0.0;

	double mu1 = 0.0, mu2 = 0.0;

	double sum_norm;

	double p_norm_F;

	ierr = this->get_G_dual(gp,E_dual_el,G_dual); CHKERRQ(ierr);

	ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);

	ierr = this->get_beta_field_at_gp(beta_el,gp,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);

	ierr = this->tranform_F_from_dual(F_0_dual,G_real,F_0_real); CHKERRQ(ierr);

	for (int k=0; k<this->dim; k++){
		for (int gamma=0; gamma<2; gamma++){ 
			F_guess[k][gamma] = F_0_real[k][gamma];
		}
	}

	// ierr = this->tranform_F_to_dual(F_guess,G_dual,F_guess_dual); CHKERRQ(ierr);

	// for (int k1=0; k1<this->dim; k1++){
	// 	for (int k2=0; k2<2; k2++){

	// 		std::cout << "F_guess is " << F_guess[k1][k2] << std::endl;
	// 	}
	// }

	while (delta_F_norm > tol_local || residual_local > tol_local)
	{	

		ierr = this->tranform_F_to_dual(F_guess,G_dual,F_guess_dual); CHKERRQ(ierr);

		ierr = this->get_C(F_guess,F_guess_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F_guess,G_real,&det_C); CHKERRQ(ierr);

		ierr = this->get_invariants(tr_C,det_C,&mu1,&mu2); CHKERRQ(ierr);

		ierr = this->get_eigenvectors_of_C(C,mu1,mu2,G_dual,eig_v1,eig_v2); CHKERRQ(ierr);

		ierr = this->first_derivative_of_invariants(eig_v1,eig_v2,G_dual,F_guess,d_mu1_F,d_mu2_F); CHKERRQ(ierr);

		ierr = this->second_derivative_of_invariants(eig_v1,eig_v2,mu1,mu2,tr_C,C,G_dual,G_real,F_guess,dd_mu1_FF,dd_mu2_FF); CHKERRQ(ierr);
		// need to code for coalescence case of eigenvalues

		sum_norm = 0.0;

		ierr = this->get_L2_norm_F_diff_F_0(F_guess, F_guess_dual, F_0_real, F_0_dual, &sum_norm); CHKERRQ(ierr);


		for (int k=0; k<this->dim; k++){

			for (int gamma=0; gamma<2; gamma++){

				index1 = k*2 + gamma;

				ierr = this->get_A_field_at_gp(A_el,index1,gp,&A_at_gp); CHKERRQ(ierr);

				// rhs_vector local
				f_local[index1] = -(1.0/this->pf)*(A_at_gp  
								+ d_mu1_F[index1]*beta_1_gp
								+ d_mu2_F[index1]*beta_2_gp
								+ (this->pf*std::pow(sum_norm,PP-1)*(F_guess_dual[k][gamma] - F_0_dual[k][gamma]))
								+ this->pf*(F_guess_dual[k][gamma] - F_0_dual[k][gamma]));


				// if (local_step < 2)
				// {
				if (std::abs(2.0*mu1 - tr_C) > 1e-3){
					
					for (int p=0; p<this->dim; p++){

						for (int omega=0;omega<2;omega++){

							index2 = p*2 + omega;

							// matrix local
							if (std::abs(sum_norm) > 1e-3){
								p_norm_F = std::pow(sum_norm,PP-2);
							}
							else {
								p_norm_F = 0.0;
							}

							k_local[index1][index2] = (1.0/this->pf)*( 
								(this->pf*std::pow(sum_norm,PP-1)*II[k][p]*G_dual[omega][gamma]) +
								(2.0*(PP-1.0)*this->pf*p_norm_F*
								(F_guess_dual[p][omega] - F_0_dual[p][omega])*(F_guess_dual[k][gamma] - F_0_dual[k][gamma]))
								+ (this->pf*II[p][k]*G_dual[omega][gamma]) +
								+ dd_mu1_FF[index1][index2]*beta_1_gp + dd_mu2_FF[index1][index2]*beta_2_gp);
						}
					}
				}
				else {

					for (int p=0; p<this->dim; p++){

						for (int omega=0;omega<2;omega++){

							index2 = p*2 + omega;

							if (std::abs(sum_norm) > 1e-3){
								p_norm_F = std::pow(sum_norm,PP-2);
							}
							else {
								p_norm_F = 0.0;
							}

							// matrix local
							k_local[index1][index2] = (1.0/this->pf)*(
								(this->pf*std::pow(sum_norm,PP-1)*II[k][p]*G_dual[omega][gamma]) + 
								(2.0*(PP-1.0)*this->pf*p_norm_F*
								(F_guess_dual[p][omega] - F_0_dual[p][omega])*(F_guess_dual[k][gamma] - F_0_dual[k][gamma]))
								+ (this->pf*II[p][k]*G_dual[omega][gamma]) +
								+ dd_mu1_FF[index1][index2]*beta_1_gp + dd_mu2_FF[index1][index2]*beta_2_gp);
						}
					}
				}
					
				// }
			}

		}

		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_local,delta_F); CHKERRQ(ierr);

		// ierr = common_utilities->petsc_local_linear_system(F_DOF,k_local,f_local,delta_F); CHKERRQ(ierr);

		// delta_F_norm and residual_local
		delta_F_norm = 0.0; 
		residual_local = 0.0;

		for (int k1=0; k1<F_DOF; k1++){
			 
			delta_F_norm = delta_F_norm +  delta_F[k1]*delta_F[k1];
			residual_local =  residual_local + f_local[k1]*f_local[k1];
		}

		delta_F_norm = std::sqrt(delta_F_norm);
		residual_local = std::sqrt(residual_local);

		// new guess
		for (int k1=0; k1<dim; k1++){
			for (int k2=0; k2<2; k2++){

				index1 = k1*2+ k2;
				F_guess[k1][k2] = F_guess[k1][k2] + delta_F[index1];
				 
			}
		}

		// for (int k1=0; k1<this->dim; k1++){
		// 	for (int k2=0; k2<2; k2++){

		// 		std::cout << "F_guess is " << F_guess[k1][k2] << std::endl;
		// 	}
		// }

		local_step = local_step + 1;

		if (local_step > 30 || residual_local > 1e+4){ //|| local_step > 500
			std::cout << "********************" << std::endl;
			std::cout << "mu1 and mu2 are " <<  mu1 << " and " << mu2 << std::endl;
			std::cout << "mu1 - mu_imp and mu2 - mu_imp are " << (mu1 - lambda_1*lambda_1) << " and " << (mu2- lambda_2*lambda_2) << std::endl;
			// std::cout << "sgn of mu - mu_imp are " << std::tanh((mu1 - lambda_1*lambda_1)/A_CONST) 
			// 							<< " and " << std::tanh((mu2- lambda_2*lambda_2)/A_CONST) <<  std::endl;
			std::cout << " local step for gp = " << gp << " is " << local_step << std::endl;
			std::cout << " local norm is " << delta_F_norm << " and " << residual_local << std::endl;
			for (int ind1=0; ind1<F_DOF; ind1++){
				// std::cout << " d_mu1_F and d_mu2_F are " << d_mu1_F[ind1] << " and " << d_mu2_F[ind1] << std::endl;
				for (int ind2=0;ind2<F_DOF;ind2++){
					// std::cout << " dd_mu1_FF and dd_mu2_FF for " << ind1 << " , " << ind2 << " are " << 
						// dd_mu1_FF[ind1][ind2] << " and " << dd_mu2_FF[ind1][ind2] << std::endl;
					std::cout << "k_local for ind1 = " << ind1 << " and ind2 = " << ind2 << " is " << k_local[ind1][ind2] << std::endl;
				}
			}
			std::cout << "********************" << std::endl;

		}
		// std::cout << "********************" << std::endl;
		// std::cout << " local step for gp = " << gp << " is " << local_step << std::endl;
		// std::cout << " local norm is " << std::fixed << std::setprecision(15) << delta_F_norm << " and " << residual_local << std::endl;
		// std::cout << "********************" << std::endl;

		// if (local_step > 1){
		// 	std::cout << "local norms are " << delta_F_norm  << " and " << residual_local << std::endl;
		// 	std::cout << "doing quasi-newton" << std::endl;
		// }

		if (local_step > 50){
			std::cout << "local norms are " << delta_F_norm  << " and " << residual_local << std::endl;
			*solution = 0;
			// std::cout << "moving to Picards solve instead of newton solve" << std::endl;
			// break;
			throw "local_step for gp convergence is greater than 50";
		}


	} // while loop ends

	// std::cout << "********************" << std::endl;
	// std::cout << " local step for gp = " << gp << " is " << local_step << std::endl;
	// std::cout << " local norm is " << delta_F_norm << " and " << residual_local << std::endl;
	// std::cout << "********************" << std::endl;

	// storing the final value to return
	if (*solution == 1){
		for (int k1=0; k1<dim; k1++){
			for (int k2=0; k2<2; k2++){

				F_sol[k1][k2] = F_guess[k1][k2];
				 
			}
		}
	}

	return 0;

}

PetscErrorCode Dual_Solve::global_newton_solve()
{

	PetscErrorCode ierr;

	PetscReal rnorm;
	PetscInt its;

	MyMat<double> beta_el(fem->node_per_el,MyVec<double> (grid->BETA_DOF));

	MyMat<double> A_el(fem->node_per_el, MyVec<double>(grid->A_DOF));

	// std::cout << " print dof_per_el. for " << fem->dof_per_el << std::endl;

	int node_el;

	double rhs_el[fem->dof_per_el];

	double **ke;

	ke = new double* [fem->dof_per_el];
	for (int i=0;i<fem->dof_per_el;i++)
		ke[i]= new double [fem->dof_per_el];

	double ke_vec[fem->dof_per_el*fem->dof_per_el];

	int indc[fem->dof_per_el];

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> x_0_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> F_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> F_guess_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> F_0_dual_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	MyMat <double> E_real_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	MyMat <double> xdef_el(fem->node_per_el,MyVec<double>(this->dim));

	MyMat <double> Eig1_el(fem->ngp, MyVec <double> (2));

	MyMat <double> Eig2_el(fem->ngp, MyVec <double> (2));

	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initializing rhs_el to be zero (passed as pointers: should not sum)
		for (int i=0;i<fem->dof_per_el;i++){
			for (int j=0;j<fem->dof_per_el;j++){
				ke[i][j] = 0.0;
			}
			rhs_el[i] = 0.0;
		}

		// getting dual fields for each elem to get mapping 
		for (int i=0;i<fem->node_per_el;i++){

			node_el = grid->el_conn[ie][i];
			
			for (int k=0;k<grid->BETA_DOF;k++){
				beta_el[i][k] = beta[node_el][k];
			}

			for (int k=0;k<grid->A_DOF;k++){
				A_el[i][k] = A_dual[node_el][k];
			}

			for (int k=0;k<this->dim;k++){
				xdef_el[i][k] = x_def[node_el][k];
			}
		}

		// getting basis vectors and jacobian for each elem
		for (int i=0;i<fem->ngp;i++){
			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
				E_real_el[i][k] = fem->E_real[ie][i][k];
				//std::cout << "E_dual_el for i and k is " << E_dual_el[i][k] << std::endl;
			}
			
			detJ_e[i] = fem->detJ_G[ie][i];	
		}

		// std::cout<< "********************" << std::endl;
		// std::cout<< "element no is " << ie << std::endl;
		// std::cout<< "********************" << std::endl;

		// getting indicies for assembly of elem vector to global vector
		for (int k=0;k<(grid->A_DOF + grid->BETA_DOF)*fem->node_per_el;k++){
			indc[k] = grid->dof_el_conn[ie][k];
		}

		// storing primal fields at gpts from previous iteration as guess and base state values (fixed)
		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<F_DOF;k++){
				F_guess_gp_el[i][k] = F_gp[ie][i][k];
			}

			for (int k=0;k<F_DOF;k++){
				F_0_dual_gp_el[i][k] = F_0_dual_gp[ie][i][k];
			}

			for (int k=0;k<this->dim;k++){
				x_0_gp_el[i][k] = x_0_gp[ie][i][k];
			}
		}

		// std::cout << " code runs up to get_primal_fields " << std::endl;

		// primal fields (obtained from mapping from dual to primal)
		ierr = this->get_primal_fields(beta_el,A_el,E_dual_el,xdef_el,F_guess_gp_el,
												x_0_gp_el,x_gp_el,F_0_dual_gp_el,F_gp_el); CHKERRQ(ierr);

		// std::cout << " code runs up to elem_stiff_force_vector_global_newton " << std::endl;

		// stiff matrix and rhs vector per elem
		ierr = this->elem_stiff_force_vector_global_newton(beta_el,A_el,E_dual_el,E_real_el,detJ_e,
											rhs_el,ke,x_gp_el,F_gp_el,F_0_dual_gp_el,Eig1_el,Eig2_el); CHKERRQ(ierr);


		// storing primal fields at gpts
		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<dim;k++){
				x_gp[ie][i][k] = x_gp_el[i][k];
			}

			for (int k=0;k<F_DOF;k++){
				F_gp[ie][i][k] = F_gp_el[i][k];
			}

			for (int k=0; k<2; k++){
				
				Eig_1_gp[ie][i][k] = Eig1_el[i][k];
				Eig_2_gp[ie][i][k] = Eig2_el[i][k];
			}
		}

		// adding values to global matrix and vector
		for (int i=0;i<fem->dof_per_el;i++){
			for (int j=0;j<fem->dof_per_el;j++){

				ke_vec[i*fem->dof_per_el+j] = ke[i][j];
			}
		}	

		// assembly to global rhs vector and stiff matrix
		ierr = VecSetValues(rhs_vec,fem->dof_per_el,indc,rhs_el,ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(K_dual,fem->dof_per_el,indc,fem->dof_per_el,indc,ke_vec,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector (due to multiple proccess)
	ierr = VecAssemblyBegin(rhs_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs_vec); CHKERRQ(ierr);

	// final assembly of global matrix (due to multiple proccess)
	ierr = MatAssemblyBegin(K_dual,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K_dual,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// boundary conditions for newton solve
	// if (this->count ==0){
	ierr = this->global_newton_solve_bcs(K_dual,rhs_vec,delta_dual_vec); CHKERRQ(ierr);
	// }

	// creating the solver context
	if (this->count == 0){
		ierr = common_utilities->ksp_mumps_solver_petsc(ksp_dual,K_dual,pc_dual,FF_dual); CHKERRQ(ierr);
	}
		
	//Solve the linear system
	ierr = KSPSolve(ksp_dual,rhs_vec,delta_dual_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(ksp_dual,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(ksp_dual,&its); CHKERRQ(ierr);

	if (::rank == 0){
		std::cout << "global newton solve residual for step no " << this->step_no << " is " << rnorm << std::endl;
	}

	// delete objects and matrices
	for (int i = 0; i<fem->dof_per_el; i++)
		delete[] ke[i];

	delete[] ke;

	return 0;
}

PetscErrorCode Dual_Solve::elem_stiff_force_vector_global_newton(MyMat <double> beta_el, MyMat <double> A_el, 
											 MyMat <double> E_dual_el, MyMat <double> E_real_el, MyVec<double> detJ_e, 
											double* rhs_el, double **ke, MyMat <double> x_gp_el,  
											MyMat <double> F_gp_el, MyMat <double> F_0_dual_gp_el,
											MyMat <double> &Eig1_el, MyMat<double> &Eig2_el)
{

	PetscErrorCode ierr;

	MyMat <double> k_local(F_DOF,MyVec <double>(F_DOF));

	// std::cout << "F_DOF is " << F_DOF << std::endl;

	MyVec <double> f_beta_1(F_DOF);
	MyVec <double> f_beta_2(F_DOF);
	MyVec <double> f_A11(F_DOF);
	MyVec <double> f_A12(F_DOF);
	MyVec <double> f_A21(F_DOF);
	MyVec <double> f_A22(F_DOF);
	MyVec <double> f_A31(F_DOF);
	MyVec <double> f_A32(F_DOF);

	MyVec <double> dF_dbeta_1(F_DOF);
	MyVec <double> dF_dbeta_2(F_DOF);
	MyVec <double> dF_dA11(F_DOF);
	MyVec <double> dF_dA12(F_DOF);
	MyVec <double> dF_dA21(F_DOF);
	MyVec <double> dF_dA22(F_DOF);
	MyVec <double> dF_dA31(F_DOF);
	MyVec <double> dF_dA32(F_DOF);

	MyMat<double> F(dim,MyVec<double>(2));

	MyMat<double> F_dual(dim,MyVec<double>(2));

	MyMat<double> F_0_dual(dim,MyVec<double>(2));

	MyMat<double> F_0_real(dim,MyVec<double>(2));

	MyMat<double> G_dual(2,MyVec<double>(2));

	MyMat<double> G_real(2,MyVec<double>(2));

	MyMat<double> C(2,MyVec<double>(2));

	double tr_C, det_C;

	double mu1, mu2;

	double sum_norm;

	double p_norm_F;

	MyVec <double> eig_v1(2);

	MyVec <double> eig_v2(2);

	MyVec <double> d_mu1_F(F_DOF);

	MyVec <double> d_mu2_F(F_DOF);

	MyMat <double> dd_mu1_FF(F_DOF,MyVec <double> (F_DOF));

	MyMat <double> dd_mu2_FF(F_DOF,MyVec <double> (F_DOF));

	double d_mu1_beta_1 = 0.0;
	double d_mu1_beta_2 = 0.0;
	double d_mu1_A11 = 0.0;
	double d_mu1_A12 = 0.0;
	double d_mu1_A21 = 0.0;
	double d_mu1_A22 = 0.0;
	double d_mu1_A31 = 0.0;
	double d_mu1_A32 = 0.0;

	double d_mu2_beta_1 = 0.0;
	double d_mu2_beta_2 = 0.0;
	double d_mu2_A11 = 0.0;
	double d_mu2_A12 = 0.0;
	double d_mu2_A21 = 0.0;
	double d_mu2_A22 = 0.0;
	double d_mu2_A31 = 0.0;
	double d_mu2_A32 = 0.0;

	int tmp, ind1, ind2;

	double beta_1_gp = 0.0, beta_2_gp = 0.0;

	double fe[fem->dof_per_el];

	double K[fem->dof_per_el][fem->dof_per_el];

	MyVec <double> beta_1_zeta(2);
	MyVec <double> beta_2_zeta(2);

	MyMat <double> dA_zeta(grid->A_DOF, MyVec <double> (2));

	MyVec <double> A_el_gp(grid->A_DOF);

	MyVec <double> Eig1(2);
	MyVec <double> Eig2(2);

	MyVec <double> A11_res1(fem->node_per_el);
	MyVec <double> A12_res1(fem->node_per_el);
	MyVec <double> A21_res1(fem->node_per_el);
	MyVec <double> A22_res1(fem->node_per_el);
	MyVec <double> A31_res1(fem->node_per_el);
	MyVec <double> A32_res1(fem->node_per_el);

	MyVec <double> A11_res2(fem->node_per_el);
	MyVec <double> A12_res2(fem->node_per_el);
	MyVec <double> A21_res2(fem->node_per_el);
	MyVec <double> A22_res2(fem->node_per_el);
	MyVec <double> A31_res2(fem->node_per_el);
	MyVec <double> A32_res2(fem->node_per_el);

	MyVec <double> A11_res3(fem->node_per_el);
	MyVec <double> A12_res3(fem->node_per_el);
	MyVec <double> A21_res3(fem->node_per_el);
	MyVec <double> A22_res3(fem->node_per_el);
	MyVec <double> A31_res3(fem->node_per_el);
	MyVec <double> A32_res3(fem->node_per_el);

	MyVec <double> A11_res4(fem->node_per_el);
	MyVec <double> A12_res4(fem->node_per_el);
	MyVec <double> A21_res4(fem->node_per_el);
	MyVec <double> A22_res4(fem->node_per_el);
	MyVec <double> A31_res4(fem->node_per_el);
	MyVec <double> A32_res4(fem->node_per_el);

	MyMat <double> dE1_zeta(grid->X_DOF, MyVec <double> (2));
	MyMat <double> dE2_zeta(grid->X_DOF, MyVec <double> (2));

	MyVec <double> A11_stiff1(grid->A_DOF);
	MyVec <double> A12_stiff1(grid->A_DOF);
	MyVec <double> A21_stiff1(grid->A_DOF);
	MyVec <double> A22_stiff1(grid->A_DOF);
	MyVec <double> A31_stiff1(grid->A_DOF);
	MyVec <double> A32_stiff1(grid->A_DOF);

	MyMat <double> A11_stiff2(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A12_stiff2(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A21_stiff2(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A22_stiff2(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A31_stiff2(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A32_stiff2(grid->A_DOF, MyVec <double> (fem->node_per_el));

	MyMat <double> A11_stiff3(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A12_stiff3(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A21_stiff3(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A22_stiff3(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A31_stiff3(grid->A_DOF, MyVec <double> (fem->node_per_el));
	MyMat <double> A32_stiff3(grid->A_DOF, MyVec <double> (fem->node_per_el));

	double dG11_zeta2 = 0.0;
	double dG12_zeta1 = 0.0;
	double dG21_zeta2 = 0.0;
	double dG22_zeta1 = 0.0;

	// derivatives of E_real wrt zeta
	for (int k=0;k<grid->X_DOF;k++){
		for (int p=0;p<fem->ngp;p++){
			
			dE1_zeta[k][0] = dE1_zeta[k][0] + fem->dpsi_c[p][0]*E_real_el[p][k];
			dE1_zeta[k][1] = dE1_zeta[k][1] + fem->dpsi_c[p][1]*E_real_el[p][k];

			dE2_zeta[k][0] = dE2_zeta[k][0] + fem->dpsi_c[p][0]*E_real_el[p][grid->X_DOF + k];
			dE2_zeta[k][1] = dE2_zeta[k][1] + fem->dpsi_c[p][1]*E_real_el[p][grid->X_DOF + k];
		}
	} 

	// summation over gauss points
	for (int i=0;i<fem->ngp;i++)
	{	

		// std::cout << " gp is " << i << std::endl;

		for (int k1=0;k1<this->dim;k1++){
			for(int k2=0;k2<2;k2++){

				tmp = k1*2 + k2;
				F[k1][k2] = F_gp_el[i][tmp];
				F_0_dual[k1][k2] = F_0_dual_gp_el[i][tmp];
			}
		}

		for (int k1=0;k1<fem->dof_per_el;k1++){
			for (int k2=0;k2<fem->dof_per_el;k2++){
				K[k1][k2] = 0.0;
			}
			fe[k1] = 0.0;
		}

		// derivative of A wrt zeta at gpts and A at gpts
		for (int k=0; k<grid->A_DOF; k++){

			dA_zeta[k][0] = 0.0;
			dA_zeta[k][1] = 0.0;
			A_el_gp[k] = 0.0;

			for (int p=0; p<fem->node_per_el; p++){
				dA_zeta[k][0] = dA_zeta[k][0] + fem->dpsi[i][p][0]*A_el[p][k];
				dA_zeta[k][1] = dA_zeta[k][1] + fem->dpsi[i][p][1]*A_el[p][k];

				A_el_gp[k] = A_el_gp[k] + fem->psi[i][p]*A_el[p][k];
			}
		}

		ierr = this->get_G_dual(i,E_dual_el,G_dual); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);		

		ierr = this->tranform_F_to_dual(F,G_dual,F_dual); CHKERRQ(ierr);

		ierr = this->tranform_F_from_dual(F_0_dual, G_real, F_0_real); CHKERRQ(ierr);

		ierr = this->get_C(F,F_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F,G_real,&det_C); CHKERRQ(ierr);
		
		ierr = this->get_invariants(tr_C,det_C,&mu1,&mu2); CHKERRQ(ierr);

		ierr = this->get_eigenvectors_of_C(C,mu1,mu2,G_dual,eig_v1,eig_v2); CHKERRQ(ierr);

		if (problem_type == inverse){
			// std::cout << "code runs upto finding eigenvectors on target shape" << std::endl;
			ierr = this->eigenvectors_on_target_shape(eig_v1,eig_v2,F,G_dual,Eig1,Eig2); CHKERRQ(ierr);
		}
		else if (problem_type == forward){
			ierr = this->eigenvectors_on_reference_shape(eig_v1,eig_v2,E_dual_el,i,Eig1,Eig2); CHKERRQ(ierr);
		}
		
		// std::cout << "code runs upto finding first_derivative_of_invariants" << std::endl;
		ierr = this->first_derivative_of_invariants(eig_v1,eig_v2,G_dual,F,d_mu1_F,d_mu2_F); CHKERRQ(ierr);

		// std::cout << "code runs upto finding second_derivative_of_invariants" << std::endl;
		ierr = this->second_derivative_of_invariants(eig_v1,eig_v2,mu1,mu2,tr_C,C,G_dual,G_real,F,dd_mu1_FF,dd_mu2_FF); CHKERRQ(ierr);

		// std::cout << "code runs upto finding get_beta_field_at_gp" << std::endl;
		ierr = this->get_beta_field_at_gp(beta_el,i,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);

		ierr = this->get_derivative_of_beta_at_gp(beta_el,i,beta_1_zeta,beta_2_zeta); CHKERRQ(ierr);

		ierr = this->get_L2_norm_F_diff_F_0(F, F_dual, F_0_real, F_0_dual, &sum_norm); CHKERRQ(ierr);

		// std::cout << "mu1 is " << mu1 << std::endl;
		// std::cout << "mu2 is " << mu2 << std::endl;

		for (int k=0; k<2; k++){

			Eig1_el[i][k] = Eig1[k];
			Eig2_el[i][k] = Eig2[k];
		}
		

		for (int k=0;k<this->dim;k++){
			for (int gamma=0;gamma<2;gamma++){

				ind1 = k*2 + gamma;

				f_beta_1[ind1] = -(1.0/this->pf)*d_mu1_F[ind1];

				f_beta_2[ind1] = -(1.0/this->pf)*d_mu2_F[ind1];

				if (ind1 == 0){
					f_A11[ind1] = -(1.0/this->pf)*1.0; f_A12[ind1] = 0.0; f_A21[ind1] = 0.0; 
					f_A22[ind1] = 0.0; 
				}
				else if (ind1 == 1){
					f_A11[ind1] = 0.0; f_A12[ind1] = -(1.0/this->pf)*1.0; f_A21[ind1] = 0.0; 
					f_A22[ind1] = 0.0; 
				}
				else if (ind1 == 2){
					f_A11[ind1] = 0.0; f_A12[ind1] = 0.0; f_A21[ind1] = -(1.0/this->pf)*1.0; f_A22[ind1] = 0.0;
				}
				else if (ind1 == 3){
					f_A11[ind1] = 0.0; f_A12[ind1] = 0.0; f_A21[ind1] = 0.0; f_A22[ind1] = -(1.0/this->pf)*1.0;
				}

				if (problem_type == forward)
				{
					if (ind1 == 0){
						f_A31[ind1] = 0.0; f_A32[ind1] = 0.0;
					}
					else if (ind1 == 1){
						f_A31[ind1] = 0.0; f_A32[ind1] = 0.0;
					}
					else if (ind1 == 2){
						f_A31[ind1] = 0.0; f_A32[ind1] = 0.0;
					}
					else if (ind1 == 3){
						f_A31[ind1] = 0.0; f_A32[ind1] = 0.0;
					}
					else if (ind1 == 4){
						f_A11[ind1] = 0.0; f_A12[ind1] = 0.0; f_A21[ind1] = 0.0; 
						f_A22[ind1] = 0.0; f_A31[ind1] = -(1.0/this->pf)*1.0; f_A32[ind1] = 0.0;
					}
					else if (ind1 == 5){
						f_A11[ind1] = 0.0; f_A12[ind1] = 0.0; f_A21[ind1] = 0.0; 
						f_A22[ind1] = 0.0; f_A31[ind1] = 0.0; f_A32[ind1] = -(1.0/this->pf)*1.0;
					}

				}
											
				if (std::abs(mu1 - mu2) > 1e-3){

					for (int p=0; p<this->dim; p++){
						for (int omega=0; omega<2; omega++){

							ind2 = p*2 + omega;

							if (std::abs(sum_norm) > 1e-3){
								p_norm_F = std::pow(sum_norm,PP-2);
							}
							else {
								p_norm_F = 0.0;
							}

							k_local[ind1][ind2] = (1.0/this->pf)*(dd_mu1_FF[ind1][ind2]*beta_1_gp
													+ dd_mu2_FF[ind1][ind2]*beta_2_gp +
													(this->pf*std::pow(sum_norm,PP-1)*II[p][k]*G_dual[gamma][omega]) + 
													(2.0*(PP-1.0)*this->pf*p_norm_F*
													(F_dual[p][omega] - F_0_dual[p][omega])*
													(F_dual[k][gamma] - F_0_dual[k][gamma])) +
													(this->pf*II[p][k]*G_dual[gamma][omega]));
						}
					} 
				}
				else {

					for (int p=0; p<this->dim; p++){
						for (int omega=0; omega<2; omega++){

							ind2 = p*2 + omega;

							if (std::abs(sum_norm) > 1e-3){
								p_norm_F = std::pow(sum_norm,PP-2);
							}
							else {
								p_norm_F = 0.0;
							}

							k_local[ind1][ind2] = (1.0/this->pf)*(dd_mu1_FF[ind1][ind2]*beta_1_gp
													+ dd_mu2_FF[ind1][ind2]*beta_2_gp +
													(this->pf*std::pow(sum_norm,PP-1)*II[p][k]*G_dual[gamma][omega]) + 
													(2.0*(PP-1.0)*this->pf*p_norm_F*
													(F_dual[p][omega] - F_0_dual[p][omega])*
													(F_dual[k][gamma] - F_0_dual[k][gamma]))+
													(this->pf*II[p][k]*G_dual[gamma][omega]));
						}
					} 
				}
				
			}
		}

		// solve for dF^k_gamma/dbeta_1
		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_beta_1,dF_dbeta_1); CHKERRQ(ierr);

		// solve for dF^k_gamma/dbeta_2
		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_beta_2,dF_dbeta_2); CHKERRQ(ierr);

		// solve for dF^k_gamma/dA11
		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_A11,dF_dA11); CHKERRQ(ierr);

		// solve for dF^k_gamma/dA12
		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_A12,dF_dA12); CHKERRQ(ierr);

		// solve for dF^k_gamma/dA21
		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_A21,dF_dA21); CHKERRQ(ierr);

		// solve for dF^k_gamma/dA22
		ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_A22,dF_dA22); CHKERRQ(ierr);

		if (problem_type == forward)
		{
			// solve for dF^k_gamma/dA31
			ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_A31,dF_dA31); CHKERRQ(ierr);

			// solve for dF^k_gamma/dA32
			ierr = common_utilities->mysolve_local_linear_system(F_DOF,k_local,f_A32,dF_dA32); CHKERRQ(ierr);
		}

		ierr = fem->dot_product(d_mu1_F,dF_dbeta_1,&d_mu1_beta_1); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu1_F,dF_dbeta_2,&d_mu1_beta_2); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu1_F,dF_dA11,&d_mu1_A11); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu1_F,dF_dA12,&d_mu1_A12); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu1_F,dF_dA21,&d_mu1_A21); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu1_F,dF_dA22,&d_mu1_A22); CHKERRQ(ierr);

		if (problem_type == forward)
		{
			ierr = fem->dot_product(d_mu1_F,dF_dA31,&d_mu1_A31); CHKERRQ(ierr);
			ierr = fem->dot_product(d_mu1_F,dF_dA32,&d_mu1_A32); CHKERRQ(ierr);
		}	


		ierr = fem->dot_product(d_mu2_F,dF_dbeta_1,&d_mu2_beta_1); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu2_F,dF_dbeta_2,&d_mu2_beta_2); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu2_F,dF_dA11,&d_mu2_A11); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu2_F,dF_dA12,&d_mu2_A12); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu2_F,dF_dA21,&d_mu2_A21); CHKERRQ(ierr);
		ierr = fem->dot_product(d_mu2_F,dF_dA22,&d_mu2_A22); CHKERRQ(ierr);

		if (problem_type == forward)
		{
			ierr = fem->dot_product(d_mu2_F,dF_dA31,&d_mu2_A31); CHKERRQ(ierr);
			ierr = fem->dot_product(d_mu2_F,dF_dA32,&d_mu2_A32); CHKERRQ(ierr);	
		}

		ierr = this->regularisation_residual_terms_calculation(i, E_real_el, A_el_gp, dA_zeta, 
															dE1_zeta, dE2_zeta, G_real, G_dual,
															A11_res1, A11_res2, A11_res3, A11_res4,
															A12_res1, A12_res2, A12_res3, A12_res4,
															A21_res1, A21_res2, A21_res3, A21_res4,
															A22_res1, A22_res2, A22_res3, A22_res4,
															A31_res1, A31_res2, A31_res3, A31_res4,
															A32_res1, A32_res2, A32_res3, A32_res4); CHKERRQ(ierr);	

		// std::cout << "code runs after finding regularisation_residual_terms_calculation" << std::endl;		

		ierr = this->regularisation_stiffness_terms_calculation(i, E_real_el, G_dual, dE1_zeta,dE2_zeta,
												A11_stiff1, A12_stiff1, A21_stiff1, A22_stiff1, A31_stiff1, A32_stiff1,
												A11_stiff2, A12_stiff2, A21_stiff2, A22_stiff2, A31_stiff2, A32_stiff2,
												A11_stiff3, A12_stiff3, A21_stiff3, A22_stiff3, A31_stiff3, A32_stiff3); CHKERRQ(ierr);

		// std::cout << "code runs after finding regularisation_stiffness_terms_calculation" << std::endl;	

		ierr = this->get_derivatives_of_G_with_zeta(i,E_real_el,dE1_zeta,dE2_zeta,
													&dG11_zeta2, &dG12_zeta1,
													&dG21_zeta2, &dG22_zeta1); CHKERRQ(ierr);	


		// std::cout << "code runs after finding get_derivatives_of_G_with_zeta" << std::endl;															

		// stiffness matrix and residual vector
		for (int p=0;p<fem->node_per_el;p++)
		{

			for (int q=0;q<fem->node_per_el;q++)
			{

				// beta_1 coupling
				K[p*grid->dof_per_node][q*grid->dof_per_node] = d_mu1_beta_1*fem->psi[i][p]*fem->psi[i][q] + 
																EPS1*(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

				K[p*grid->dof_per_node][q*grid->dof_per_node+1] = d_mu1_beta_2*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node][q*grid->dof_per_node+2] = d_mu1_A11*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node][q*grid->dof_per_node+3] = d_mu1_A12*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node][q*grid->dof_per_node+4] = d_mu1_A21*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node][q*grid->dof_per_node+5] = d_mu1_A22*fem->psi[i][p]*fem->psi[i][q];

				if (problem_type == forward)
				{
					K[p*grid->dof_per_node][q*grid->dof_per_node+6] = d_mu1_A31*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node][q*grid->dof_per_node+7] = d_mu1_A32*fem->psi[i][p]*fem->psi[i][q];
				}

				// ************************** //
				// beta_2 coupling
				K[p*grid->dof_per_node+1][q*grid->dof_per_node] = d_mu2_beta_1*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node+1][q*grid->dof_per_node+1] = d_mu2_beta_2*fem->psi[i][p]*fem->psi[i][q] +
																EPS2*(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

				K[p*grid->dof_per_node+1][q*grid->dof_per_node+2] = d_mu2_A11*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node+1][q*grid->dof_per_node+3] = d_mu2_A12*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node+1][q*grid->dof_per_node+4] = d_mu2_A21*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node+1][q*grid->dof_per_node+5] = d_mu2_A22*fem->psi[i][p]*fem->psi[i][q];

				if (problem_type == forward)
				{
					K[p*grid->dof_per_node+1][q*grid->dof_per_node+6] = d_mu2_A31*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+1][q*grid->dof_per_node+7] = d_mu2_A32*fem->psi[i][p]*fem->psi[i][q];
				}

				// ************************** //
				// A11 coupling
				K[p*grid->dof_per_node+2][q*grid->dof_per_node] = dF_dbeta_1[0]*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node+2][q*grid->dof_per_node+1] = dF_dbeta_2[0]*fem->psi[i][p]*fem->psi[i][q];

				K[p*grid->dof_per_node+2][q*grid->dof_per_node+2] = dF_dA11[0]*fem->psi[i][p]*fem->psi[i][q] + 
																		(-1.0/this->px)*fem->dpsi[i][p][0]*fem->dpsi[i][q][0] + 
																	EPS3*G_real[0][0]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A11_stiff1[0]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A11_stiff2[0][q])
																	+ EPS3*(fem->psi[i][q]*A11_stiff3[0][p])
																	+ ONE_BY_EPS*fem->dpsi[i][p][1]*fem->dpsi[i][q][1]
																	+ PEN_CURL*(G_real[0][0]*fem->dpsi[i][p][1] -
																		G_real[0][1]*fem->dpsi[i][p][0] + 
																		fem->psi[i][p]*dG11_zeta2 - 
																		fem->psi[i][p]*dG12_zeta1)
																		*(G_real[0][0]*fem->dpsi[i][q][1] + 
																		fem->psi[i][q]*dG11_zeta2 - 
																		G_real[0][1]*fem->dpsi[i][q][0] - 
																		fem->psi[i][q]*dG12_zeta1);
				
				K[p*grid->dof_per_node+2][q*grid->dof_per_node+3] = (-1.0/this->px)*fem->dpsi[i][p][0]*fem->dpsi[i][q][1]
																	+ dF_dA12[0]*fem->psi[i][p]*fem->psi[i][q] +
																	EPS3*G_real[1][0]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A11_stiff1[1]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A11_stiff2[1][q])
																	+ EPS3*(fem->psi[i][q]*A11_stiff3[1][p])
																	- ONE_BY_EPS*fem->dpsi[i][p][1]*fem->dpsi[i][q][0]
																	+ PEN_CURL*(G_real[0][0]*fem->dpsi[i][p][1] -
																		G_real[0][1]*fem->dpsi[i][p][0] + 
																		fem->psi[i][p]*dG11_zeta2 - 
																		fem->psi[i][p]*dG12_zeta1)
																		*(fem->dpsi[i][q][1]*G_real[1][0] 
																		- fem->dpsi[i][q][0]*G_real[1][1]
																		+ fem->psi[i][q]*dG21_zeta2 
																		+ fem->psi[i][q]*dG22_zeta1);
																		
				
				K[p*grid->dof_per_node+2][q*grid->dof_per_node+4] = dF_dA21[0]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+2][q*grid->dof_per_node+5] = dF_dA22[0]*fem->psi[i][p]*fem->psi[i][q];

				if (problem_type == forward)
				{
					K[p*grid->dof_per_node+2][q*grid->dof_per_node+6] = dF_dA31[0]*fem->psi[i][p]*fem->psi[i][q];
				
					K[p*grid->dof_per_node+2][q*grid->dof_per_node+7] = dF_dA32[0]*fem->psi[i][p]*fem->psi[i][q];
				}

				// ************************** //
				// A12 coupling
				K[p*grid->dof_per_node+3][q*grid->dof_per_node] = dF_dbeta_1[1]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+3][q*grid->dof_per_node+1] = dF_dbeta_2[1]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+3][q*grid->dof_per_node+2] = (-1.0/this->px)*fem->dpsi[i][p][1]*fem->dpsi[i][q][0]
																	+ dF_dA11[1]*fem->psi[i][p]*fem->psi[i][q] + 
																	EPS3*G_real[0][1]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A12_stiff1[0]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A12_stiff2[0][q])
																	+ EPS3*(fem->psi[i][q]*A12_stiff3[0][p])
																	- ONE_BY_EPS*fem->dpsi[i][p][0]*fem->dpsi[i][q][1]
																	+ PEN_CURL*(fem->dpsi[i][p][1]*G_real[1][0]
																	 + fem->psi[i][p]*dG21_zeta2 
																	 - fem->dpsi[i][p][0]*G_real[1][1]
																	 - fem->psi[i][p]*dG22_zeta1)*
																	 (fem->dpsi[i][q][1]*G_real[0][0] + 
																	 fem->psi[i][q]*dG11_zeta2 - 
																	 fem->dpsi[i][q][0]*G_real[0][1] -
																	 fem->psi[i][q]*dG12_zeta1);
																	
				
				K[p*grid->dof_per_node+3][q*grid->dof_per_node+3] = dF_dA12[1]*fem->psi[i][p]*fem->psi[i][q] +
																	(-1.0/this->px)*fem->dpsi[i][p][1]*fem->dpsi[i][q][1] +
																	EPS3*G_real[1][1]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A12_stiff1[1]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A12_stiff2[1][q])
																	+ EPS3*(fem->psi[i][q]*A12_stiff3[1][p])
																	+ ONE_BY_EPS*fem->dpsi[i][p][0]*fem->dpsi[i][q][0]
																	+ PEN_CURL*(fem->dpsi[i][p][1]*G_real[1][0]
																	 + fem->psi[i][p]*dG21_zeta2 
																	 - fem->dpsi[i][p][0]*G_real[1][1]
																	 - fem->psi[i][p]*dG22_zeta1)*
																	 (fem->dpsi[i][q][1]*G_real[1][0] + 
																	 fem->psi[i][q]*dG21_zeta2 - 
																	 fem->dpsi[i][q][0]*G_real[1][1] -
																	 fem->psi[i][q]*dG22_zeta1);
				
				K[p*grid->dof_per_node+3][q*grid->dof_per_node+4] = dF_dA21[1]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+3][q*grid->dof_per_node+5] = dF_dA22[1]*fem->psi[i][p]*fem->psi[i][q];

				if (problem_type == forward)
				{
					K[p*grid->dof_per_node+3][q*grid->dof_per_node+6] = dF_dA31[1]*fem->psi[i][p]*fem->psi[i][q];
				
					K[p*grid->dof_per_node+3][q*grid->dof_per_node+7] = dF_dA32[1]*fem->psi[i][p]*fem->psi[i][q];
				}

				// ************************** //
				// A21 coupling
				K[p*grid->dof_per_node+4][q*grid->dof_per_node] = dF_dbeta_1[2]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+4][q*grid->dof_per_node+1] = dF_dbeta_2[2]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+4][q*grid->dof_per_node+2] = dF_dA11[2]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+4][q*grid->dof_per_node+3] = dF_dA12[2]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+4][q*grid->dof_per_node+4] = dF_dA21[2]*fem->psi[i][p]*fem->psi[i][q] +
																		(-1.0/this->px)*fem->dpsi[i][p][0]*fem->dpsi[i][q][0] +
																	EPS3*G_real[0][0]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A21_stiff1[2]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A21_stiff2[2][q])
																	+ EPS3*(fem->psi[i][q]*A21_stiff3[2][p])
																	+ ONE_BY_EPS*fem->dpsi[i][p][1]*fem->dpsi[i][q][1]
																	+ PEN_CURL*(fem->dpsi[i][p][1]*G_real[0][0]
																		+ fem->psi[i][p]*dG11_zeta2
																		- fem->dpsi[i][p][0]*G_real[0][1]
																		+ fem->psi[i][p]*dG12_zeta1)
																		*(fem->dpsi[i][q][1]*G_real[0][0]
																		+ fem->psi[i][q]*dG11_zeta2
																		- fem->dpsi[i][q][0]*G_real[0][1]
																		- fem->psi[i][q]*dG12_zeta1);
				
				K[p*grid->dof_per_node+4][q*grid->dof_per_node+5] =  (-1.0/this->px)*fem->dpsi[i][p][0]*fem->dpsi[i][q][1] 
																	+ dF_dA22[2]*fem->psi[i][p]*fem->psi[i][q] + 
																	EPS3*G_real[1][0]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A21_stiff1[3]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A21_stiff2[3][q])
																	+ EPS3*(fem->psi[i][q]*A21_stiff3[3][p])
																	- ONE_BY_EPS*fem->dpsi[i][p][1]*fem->dpsi[i][q][0]
																	+ PEN_CURL*(fem->dpsi[i][p][1]*G_real[0][0]
																		+ fem->psi[i][p]*dG11_zeta2
																		- fem->dpsi[i][p][0]*G_real[0][1]
																		+ fem->psi[i][p]*dG12_zeta1)
																		*(fem->dpsi[i][q][1]*G_real[1][0]
																		+ fem->psi[i][q]*dG21_zeta2
																		- fem->dpsi[i][q][0]*G_real[1][1]
																		- fem->psi[i][q]*dG22_zeta1);
				
				if (problem_type == forward)
				{
					K[p*grid->dof_per_node+4][q*grid->dof_per_node+6] = dF_dA31[2]*fem->psi[i][p]*fem->psi[i][q];
				
					K[p*grid->dof_per_node+4][q*grid->dof_per_node+7] = dF_dA32[2]*fem->psi[i][p]*fem->psi[i][q];
				}
																		
				// ************************** //
				// A22 coupling
				K[p*grid->dof_per_node+5][q*grid->dof_per_node] = dF_dbeta_1[3]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+5][q*grid->dof_per_node+1] = dF_dbeta_2[3]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+5][q*grid->dof_per_node+2] = dF_dA11[3]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+5][q*grid->dof_per_node+3] = dF_dA12[3]*fem->psi[i][p]*fem->psi[i][q];
				
				K[p*grid->dof_per_node+5][q*grid->dof_per_node+4] = (-1.0/this->px)*fem->dpsi[i][p][1]*fem->dpsi[i][q][0]
																	+ dF_dA21[3]*fem->psi[i][p]*fem->psi[i][q] + 
																	EPS3*G_real[0][1]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A22_stiff1[2]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A22_stiff2[2][q])
																	+ EPS3*(fem->psi[i][q]*A22_stiff3[2][p])
																	- ONE_BY_EPS*fem->dpsi[i][p][0]*fem->dpsi[i][q][1]
																	+ PEN_CURL*(fem->dpsi[i][p][1]*G_real[1][0]
																	  + fem->psi[i][p]*dG21_zeta2
																	  - fem->dpsi[i][p][0]*G_real[1][1]
																	  - fem->psi[i][p]*dG22_zeta1)
																	  *(fem->dpsi[i][q][1]*G_real[0][0]
																	  + fem->psi[i][q]*dG11_zeta2
																	  - fem->dpsi[i][q][0]*G_real[0][1]
																	  - fem->psi[i][q]*dG12_zeta1);
																		
				
				K[p*grid->dof_per_node+5][q*grid->dof_per_node+5] = dF_dA22[3]*fem->psi[i][p]*fem->psi[i][q] +
																	(-1.0/this->px)*fem->dpsi[i][p][1]*fem->dpsi[i][q][1] + 
																	EPS3*G_real[1][1]*
																	(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
																	fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] + 
																	fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																	
																	+ EPS3*(A22_stiff1[3]*fem->psi[i][p]*fem->psi[i][q])
																	+ EPS3*(fem->psi[i][p]*A22_stiff2[3][q])
																	+ EPS3*(fem->psi[i][q]*A22_stiff3[3][p])
																	+ ONE_BY_EPS*fem->dpsi[i][p][0]*fem->dpsi[i][q][0]
																	+ PEN_CURL*(fem->dpsi[i][p][1]*G_real[1][0]
																	  + fem->psi[i][p]*dG21_zeta2
																	  - fem->dpsi[i][p][0]*G_real[1][1]
																	  - fem->psi[i][p]*dG22_zeta1)
																	  *(fem->dpsi[i][q][1]*G_real[1][0]
																	  + fem->psi[i][q]*dG21_zeta2
																	  - fem->dpsi[i][q][0]*G_real[1][1]
																	  - fem->psi[i][q]*dG22_zeta1);

				if (problem_type == forward)
				{
					K[p*grid->dof_per_node+5][q*grid->dof_per_node+6] = dF_dA31[3]*fem->psi[i][p]*fem->psi[i][q];
				
					K[p*grid->dof_per_node+5][q*grid->dof_per_node+7] = dF_dA32[3]*fem->psi[i][p]*fem->psi[i][q];
				}
				// ************************** //

				if (problem_type == forward)
				{	
					// ************************** //
					// A31 coupling
					K[p*grid->dof_per_node+6][q*grid->dof_per_node] = dF_dbeta_1[4]*fem->psi[i][p]*fem->psi[i][q];
					
					K[p*grid->dof_per_node+6][q*grid->dof_per_node+1] = dF_dbeta_2[4]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+6][q*grid->dof_per_node+2] = dF_dA11[4]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+6][q*grid->dof_per_node+3] = dF_dA12[4]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+6][q*grid->dof_per_node+4] = dF_dA21[4]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+6][q*grid->dof_per_node+5] = dF_dA22[4]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+6][q*grid->dof_per_node+6] = dF_dA31[4]*fem->psi[i][p]*fem->psi[i][q] + 
																		(-1.0/this->px)*fem->dpsi[i][p][0]*fem->dpsi[i][q][0] +
																		EPS3*G_real[0][0]*
																		(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] +
																		fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
																		fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
																		fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																		+ EPS3*(A31_stiff1[4]*fem->psi[i][p]*fem->psi[i][q])
																		+ EPS3*(fem->psi[i][p]*A31_stiff2[4][q])
																		+ EPS3*(fem->psi[i][q]*A31_stiff3[4][p]);

					K[p*grid->dof_per_node+6][q*grid->dof_per_node+7] = dF_dA32[4]*fem->psi[i][p]*fem->psi[i][q] + 
																		(-1.0/this->px)*fem->dpsi[i][p][0]*fem->dpsi[i][q][1] +
																		EPS3*G_real[1][0]*
																		(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] +
																		fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
																		fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
																		fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																		+ EPS3*(A31_stiff1[5]*fem->psi[i][p]*fem->psi[i][q])
																		+ EPS3*(fem->psi[i][p]*A31_stiff2[5][q])
																		+ EPS3*(fem->psi[i][q]*A31_stiff3[5][p]);

					// ************************** //
					// A32 coupling
					K[p*grid->dof_per_node+7][q*grid->dof_per_node] = dF_dbeta_1[5]*fem->psi[i][p]*fem->psi[i][q];
					
					K[p*grid->dof_per_node+7][q*grid->dof_per_node+1] = dF_dbeta_2[5]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+7][q*grid->dof_per_node+2] = dF_dA11[5]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+7][q*grid->dof_per_node+3] = dF_dA12[5]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+7][q*grid->dof_per_node+4] = dF_dA21[5]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+7][q*grid->dof_per_node+5] = dF_dA22[5]*fem->psi[i][p]*fem->psi[i][q];

					K[p*grid->dof_per_node+7][q*grid->dof_per_node+6] = dF_dA31[5]*fem->psi[i][p]*fem->psi[i][q] + 
																			(-1.0/this->px)*fem->dpsi[i][p][1]*fem->dpsi[i][q][0] + 
																			EPS3*G_real[0][1]*
																			(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] +
																			fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
																			fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
																			fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																			+ EPS3*(A32_stiff1[4]*fem->psi[i][p]*fem->psi[i][q])
																			+ EPS3*(fem->psi[i][p]*A32_stiff2[4][q])
																			+ EPS3*(fem->psi[i][q]*A32_stiff3[4][p]);

					K[p*grid->dof_per_node+7][q*grid->dof_per_node+7] = dF_dA32[5]*fem->psi[i][p]*fem->psi[i][q] + 
																			(-1.0/this->px)*fem->dpsi[i][p][1]*fem->dpsi[i][q][1] +
																			EPS3*G_real[1][1]*
																			(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] +
																			fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
																			fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
																			fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1])
																			+ EPS3*(A32_stiff1[5]*fem->psi[i][p]*fem->psi[i][q])
																			+ EPS3*(fem->psi[i][p]*A32_stiff2[5][q])
																			+ EPS3*(fem->psi[i][q]*A32_stiff3[5][p]);
				}

			}

			// residual vector
			//beta1
			fe[p*grid->dof_per_node] = fem->psi[i][p]*(mu1- (lambda_1*lambda_1)) 
										+ EPS1*(beta_1_zeta[0]*fem->dpsi[i][p][0]*G_dual[0][0] + 
												beta_1_zeta[1]*fem->dpsi[i][p][0]*G_dual[1][0] + 
												beta_1_zeta[0]*fem->dpsi[i][p][1]*G_dual[0][1] + 
												beta_1_zeta[1]*fem->dpsi[i][p][1]*G_dual[1][1]);

			//beta2
			fe[p*grid->dof_per_node+1] = fem->psi[i][p]*(mu2- (lambda_2*lambda_2))
											+ EPS2*(beta_2_zeta[0]*fem->dpsi[i][p][0]*G_dual[0][0] + 
												beta_2_zeta[1]*fem->dpsi[i][p][0]*G_dual[1][0] + 
												beta_2_zeta[0]*fem->dpsi[i][p][1]*G_dual[0][1] + 
												beta_2_zeta[1]*fem->dpsi[i][p][1]*G_dual[1][1]);

			// for (int k1=0;k1<dim;k1++){
			// 	for (int k2=0;k2<2;k2++){

			// 		fe[p*grid->dof_per_node+2+k1*this->dim+k2] = fem->psi[i][p]*F[k1][k2] 
			// 											+ fem->dpsi[i][p][k2]*x_gp_el[i][k1]; 
			// 	}
			// }

			//A11
			fe[p*grid->dof_per_node + 2] = fem->psi[i][p]*F[0][0] + fem->dpsi[i][p][0]*x_gp_el[i][0] 

											+ PEN_CURL*(((dA_zeta[0][1]*G_real[0][0] + dA_zeta[1][1]*G_real[1][0]
											+ A_el_gp[0]*dG11_zeta2 + A_el_gp[1]*dG21_zeta2) - 
											(dA_zeta[0][0]*G_real[0][1] + dA_zeta[1][0]*G_real[1][1]
											 + A_el_gp[0]*dG12_zeta1 + A_el_gp[1]*dG22_zeta1))*
											 (G_real[0][0]*fem->dpsi[i][p][1] - G_real[0][1]*fem->dpsi[i][p][0]
											 + dG11_zeta2*fem->psi[i][p] - dG12_zeta1*fem->psi[i][p]))
											+ ONE_BY_EPS*fem->dpsi[i][p][1]*(dA_zeta[0][1] - dA_zeta[1][0])
											
											+ EPS3_res*(A11_res1[p] + A11_res2[p] + A11_res3[p] + A11_res4[p]);
			
			// std::cout << " fe at " << p << " dof for " << i << " gp is " << fe[p*grid->dof_per_node + 2] << std::endl;
			// std::cout << " fe at " << p << " dof for " << i << " gp is " << fe[p*grid->dof_per_node + 2] << std::endl;

			//A12
			fe[p*grid->dof_per_node + 3] = fem->psi[i][p]*F[0][1] + fem->dpsi[i][p][1]*x_gp_el[i][0]
											
											+ PEN_CURL*(((dA_zeta[0][1]*G_real[0][0] + dA_zeta[1][1]*G_real[1][0]
											+ A_el_gp[0]*dG11_zeta2 + A_el_gp[1]*dG21_zeta2) - 
											(dA_zeta[0][0]*G_real[0][1] + dA_zeta[1][0]*G_real[1][1]
											 + A_el_gp[0]*dG12_zeta1 + A_el_gp[1]*dG22_zeta1))*
											 (G_real[1][0]*fem->dpsi[i][p][1] - G_real[1][1]*fem->dpsi[i][p][0]
											 + dG21_zeta2*fem->psi[i][p] - dG22_zeta1*fem->psi[i][p]))
											 - ONE_BY_EPS*fem->dpsi[i][p][0]*(dA_zeta[0][1] - dA_zeta[1][0])
											 
											 + EPS3_res*(A12_res1[p] + A12_res2[p] + A12_res3[p] +  A12_res4[p]);
			//A21
			fe[p*grid->dof_per_node + 4] = fem->psi[i][p]*F[1][0] + fem->dpsi[i][p][0]*x_gp_el[i][1]
											
											+ PEN_CURL*(((dA_zeta[2][1]*G_real[0][0] + dA_zeta[3][1]*G_real[1][0]
											+ A_el_gp[2]*dG11_zeta2 + A_el_gp[3]*dG21_zeta2) - 
											(dA_zeta[2][0]*G_real[0][1] + dA_zeta[3][0]*G_real[1][1]
											 + A_el_gp[2]*dG12_zeta1 + A_el_gp[3]*dG22_zeta1))*
											 (G_real[0][0]*fem->dpsi[i][p][1] - G_real[0][1]*fem->dpsi[i][p][0]
											 + dG11_zeta2*fem->psi[i][p] - dG12_zeta1*fem->psi[i][p]))
											+ ONE_BY_EPS*fem->dpsi[i][p][1]*(dA_zeta[2][1] - dA_zeta[3][0])
											 
											+ EPS3_res*(A21_res1[p] + A21_res2[p] + A21_res3[p] + A21_res4[p]);
			//A22
			fe[p*grid->dof_per_node + 5] = fem->psi[i][p]*F[1][1] + fem->dpsi[i][p][1]*x_gp_el[i][1]
											
											+ PEN_CURL*(((dA_zeta[2][1]*G_real[0][0] + dA_zeta[3][1]*G_real[1][0]
											+ A_el_gp[2]*dG11_zeta2 + A_el_gp[3]*dG21_zeta2) - 
											(dA_zeta[2][0]*G_real[0][1] + dA_zeta[3][0]*G_real[1][1]
											 + A_el_gp[2]*dG12_zeta1 + A_el_gp[3]*dG22_zeta1))*
											 (G_real[1][0]*fem->dpsi[i][p][1] - G_real[1][1]*fem->dpsi[i][p][0]
											 + dG21_zeta2*fem->psi[i][p] - dG22_zeta1*fem->psi[i][p]))
											- ONE_BY_EPS*fem->dpsi[i][p][0]*(dA_zeta[2][1] - dA_zeta[3][0])
										   	
											+ EPS3_res*(A22_res1[p] + A22_res2[p] + A22_res3[p] + A22_res4[p]);
			
			if (problem_type == forward){
				
				//A31
				fe[p*grid->dof_per_node + 6] = fem->psi[i][p]*F[2][0] + fem->dpsi[i][p][0]*x_gp_el[i][2]
												+ EPS3_res*(A31_res1[p] + A31_res2[p] + A31_res3[p] + A31_res4[p]);

				//A32
				fe[p*grid->dof_per_node + 7] = fem->psi[i][p]*F[2][1] + fem->dpsi[i][p][1]*x_gp_el[i][2]
												+ EPS3_res*(A32_res1[p] + A32_res2[p] + A32_res3[p] + A32_res4[p]);
			}
		}

		for (int p=0;p<fem->dof_per_el;p++){
			for (int q=0;q<fem->dof_per_el;q++){

				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[i]*fem->wt[i];
			}

			rhs_el[p] = rhs_el[p] + (-1.0)*fe[p]*detJ_e[i]*fem->wt[i];
			
			// std::cout << " fe at " << p << " dof for " << i << " gp is " << fe[p] << std::endl;
		}

	}// loop over gpts

	// std::cout << " rhs at 0 node" << " for A11 dof is " << rhs_el[0*grid->dof_per_node+2] << std::endl;
	// std::cout << " rhs at 1 node" << " for A11 dof is " << rhs_el[1*grid->dof_per_node+2] << std::endl;
	// std::cout << " rhs at 2 node" << " for A11 dof is " << rhs_el[2*grid->dof_per_node+2] << std::endl;
	// std::cout << " rhs at 3 node" << " for A11 dof is " << rhs_el[3*grid->dof_per_node+2] << std::endl;

	return 0;
}

PetscErrorCode Dual_Solve::get_derivatives_of_G_with_zeta(int i, MyMat <double> E_real_el,
													MyMat <double> dE1_zeta, MyMat <double> dE2_zeta,
													double *dG11_zeta2, double *dG12_zeta1,
													double *dG21_zeta2, double *dG22_zeta1)
{
	PetscErrorCode ierr;

	*dG11_zeta2 = 0.0;
	*dG12_zeta1 = 0.0;
	*dG21_zeta2 = 0.0;
	*dG22_zeta1 = 0.0;

	MyVec <double> v1(grid->X_DOF);
	MyVec <double> v2(grid->X_DOF);
	MyVec <double> v3(grid->X_DOF);
	MyVec <double> v4(grid->X_DOF);

	double scal1, scal2;

	for (int k=0;k<grid->X_DOF; k++)
	{
		v1[k] = dE1_zeta[k][1];
		v2[k] = E_real_el[i][k];
	}

	scal1  = 0.0;

	ierr = fem->dot_product(v1, v2, &scal1);  CHKERRQ(ierr);

	*dG11_zeta2 = 2.0*scal1;

	for (int k=0; k<grid->X_DOF; k++){
		v1[k] = dE1_zeta[k][0];
		v2[k] = E_real_el[i][grid->X_DOF + k];

		v3[k] = dE2_zeta[k][0];
		v4[k] = E_real_el[i][k];

	}

	scal1  = 0.0; scal2 = 0.0;

	ierr = fem->dot_product(v1, v2, &scal1);  CHKERRQ(ierr);

	ierr = fem->dot_product(v3, v4, &scal2);  CHKERRQ(ierr);

	*dG12_zeta1 = scal1 +  scal2;

	for (int k=0; k<grid->X_DOF; k++){
		v1[k] = dE2_zeta[k][1];
		v2[k] = E_real_el[i][k];

		v3[k] = dE1_zeta[k][1];
		v4[k] = E_real_el[i][grid->X_DOF + k];

	}

	scal1  = 0.0; scal2 = 0.0;

	ierr = fem->dot_product(v1, v2, &scal1);  CHKERRQ(ierr);

	ierr = fem->dot_product(v3, v4, &scal2);  CHKERRQ(ierr);

	*dG21_zeta2 = scal1 +  scal2;

	for (int k=0; k<grid->X_DOF; k++){
		v1[k] = dE2_zeta[k][0];
		v2[k] = E_real_el[i][grid->X_DOF + k];

	}

	scal1  = 0.0; 

	ierr = fem->dot_product(v1, v2, &scal1);  CHKERRQ(ierr);

	*dG22_zeta1 = 2.0*scal1;


	return 0;
}

PetscErrorCode Dual_Solve::regularisation_residual_terms_calculation(int i, MyMat <double> E_real_el, 
													MyVec <double> A_el_gp, MyMat <double> dA_zeta, 
													MyMat <double> dE1_zeta, MyMat <double> dE2_zeta,
													MyMat <double> G_real, MyMat <double> G_dual,
													MyVec <double> &A11_res1, MyVec <double> &A11_res2,
													MyVec <double> &A11_res3, MyVec <double> &A11_res4,
													MyVec <double> &A12_res1, MyVec <double> &A12_res2,
													MyVec <double> &A12_res3, MyVec <double> &A12_res4,
													MyVec <double> &A21_res1, MyVec <double> &A21_res2,
													MyVec <double> &A21_res3, MyVec <double> &A21_res4,
													MyVec <double> &A22_res1, MyVec <double> &A22_res2,
													MyVec <double> &A22_res3, MyVec <double> &A22_res4,
													MyVec <double> &A31_res1, MyVec <double> &A31_res2,
													MyVec <double> &A31_res3, MyVec <double> &A31_res4,
													MyVec <double> &A32_res1, MyVec <double> &A32_res2,
													MyVec <double> &A32_res3, MyVec <double> &A32_res4)
{

	PetscErrorCode ierr;

	MyVec <double> v1(grid->X_DOF);
	MyVec <double> v2(grid->X_DOF);
	MyVec <double> v3(grid->X_DOF);
	MyVec <double> v4(grid->X_DOF);

	MyVec <double> Ereal1(grid->X_DOF);
	MyVec <double> Ereal2(grid->X_DOF);

	double scal1, scal2, scal3, scal4, scal5, scal6;

	// E real
	Ereal1[0] = E_real_el[i][0];
	Ereal1[1] = E_real_el[i][1];
	Ereal1[2] = E_real_el[i][2];

	Ereal2[0] = E_real_el[i][grid->X_DOF];
	Ereal2[1] = E_real_el[i][grid->X_DOF + 1];
	Ereal2[2] = E_real_el[i][grid->X_DOF + 2];

	// residual terms calculation
	for (int p=0;p<fem->node_per_el; p++)
	{
		
		A11_res1[p] = 0.0; A12_res1[p] = 0.0; A21_res1[p] = 0.0; A22_res1[p] = 0.0;
		A11_res2[p] = 0.0; A12_res2[p] = 0.0; A21_res2[p] = 0.0; A22_res2[p] = 0.0;
		A11_res3[p] = 0.0; A12_res3[p] = 0.0; A21_res3[p] = 0.0; A22_res3[p] = 0.0;
		A11_res4[p] = 0.0; A12_res4[p] = 0.0; A21_res4[p] = 0.0; A22_res4[p] = 0.0;

		A31_res1[p] = 0.0; A32_res1[p] = 0.0; 
		A31_res2[p] = 0.0; A32_res2[p] = 0.0; 
		A31_res3[p] = 0.0; A32_res3[p] = 0.0; 
		A31_res4[p] = 0.0; A32_res4[p] = 0.0;

		for (int gamma=0; gamma<2; gamma++){
		
			for (int alpha=0; alpha<2; alpha++){
				for (int beta=0; beta<2; beta++){

					// term 1
					A11_res1[p] = A11_res1[p] + fem->dpsi[i][p][beta]*dA_zeta[gamma][alpha]*
															G_real[gamma][0]*G_dual[alpha][beta];

					A12_res1[p] = A12_res1[p] + fem->dpsi[i][p][beta]*dA_zeta[gamma][alpha]*
															G_real[gamma][1]*G_dual[alpha][beta];
					
					A21_res1[p] = A21_res1[p] + fem->dpsi[i][p][beta]*dA_zeta[2 + gamma][alpha]*
															G_real[gamma][0]*G_dual[alpha][beta];

					A22_res1[p] = A22_res1[p] + fem->dpsi[i][p][beta]*dA_zeta[2 + gamma][alpha]*
															G_real[gamma][1]*G_dual[alpha][beta];

					if (problem_type == forward)
					{
						A31_res1[p] = A31_res1[p] +  fem->dpsi[i][p][beta]*dA_zeta[4 + gamma][alpha]*
																G_real[gamma][0]*G_dual[alpha][beta];

						A32_res1[p] = A32_res1[p] +  fem->dpsi[i][p][beta]*dA_zeta[4 + gamma][alpha]*
																G_real[gamma][1]*G_dual[alpha][beta];
					}															

					// term 2
					if (gamma == 0){
						
						v1[0] = dE1_zeta[0][alpha];
						v1[1] = dE1_zeta[1][alpha];
						v1[2] = dE1_zeta[2][alpha];
					}
					else if (gamma == 1){

						v1[0] = dE2_zeta[0][alpha];
						v1[1] = dE2_zeta[1][alpha];
						v1[2] = dE2_zeta[2][alpha];
					}

					v2[0] = dE1_zeta[0][beta];
					v2[1] = dE1_zeta[1][beta];
					v2[2] = dE1_zeta[2][beta];

					v3[0] = dE2_zeta[0][beta];
					v3[1] = dE2_zeta[1][beta];
					v3[2] = dE2_zeta[2][beta];

					scal1 = 0.0; scal2 = 0.0; scal3 = 0.0;
					scal4 = 0.0; scal5 = 0.0; scal6 = 0.0;

					ierr = fem->dot_product(v1,v2,&scal1); CHKERRQ(ierr);
					ierr = fem->dot_product(v1,v3,&scal2); CHKERRQ(ierr);

					A11_res2[p] = A11_res2[p] + fem->psi[i][p]*A_el_gp[gamma]*G_dual[alpha][beta]*scal1;
					A12_res2[p] = A12_res2[p] + fem->psi[i][p]*A_el_gp[gamma]*G_dual[alpha][beta]*scal2;
					A21_res2[p] = A21_res2[p] + fem->psi[i][p]*A_el_gp[2 + gamma]*G_dual[alpha][beta]*scal1;
					A22_res2[p] = A22_res2[p] + fem->psi[i][p]*A_el_gp[2 + gamma]*G_dual[alpha][beta]*scal2;
					
					if (problem_type == forward)
					{
						A31_res2[p] = A31_res2[p] + fem->psi[i][p]*A_el_gp[4 + gamma]*G_dual[alpha][beta]*scal1;
						A32_res2[p] = A32_res2[p] + fem->psi[i][p]*A_el_gp[4 + gamma]*G_dual[alpha][beta]*scal2;
					}
					
					v4[0] = E_real_el[i][grid->X_DOF*gamma];
					v4[1] = E_real_el[i][grid->X_DOF*gamma+1];
					v4[2] = E_real_el[i][grid->X_DOF*gamma+2];

					ierr = fem->dot_product(v4,v2,&scal3); CHKERRQ(ierr);
					ierr = fem->dot_product(v4,v3,&scal4); CHKERRQ(ierr);

					//term 3
					A11_res3[p] = A11_res3[p] + fem->psi[i][p]*dA_zeta[gamma][alpha]*G_dual[alpha][beta]*scal3;
					A12_res3[p] = A12_res3[p] + fem->psi[i][p]*dA_zeta[gamma][alpha]*G_dual[alpha][beta]*scal4;
					A21_res3[p] = A21_res3[p] + fem->psi[i][p]*dA_zeta[2 + gamma][alpha]*G_dual[alpha][beta]*scal3;
					A22_res3[p] = A22_res3[p] + fem->psi[i][p]*dA_zeta[2 + gamma][alpha]*G_dual[alpha][beta]*scal4;

					if (problem_type == forward)
					{
						A31_res3[p] = A31_res3[p] + fem->psi[i][p]*dA_zeta[4 + gamma][alpha]*G_dual[alpha][beta]*scal3;
						A32_res3[p] = A32_res3[p] + fem->psi[i][p]*dA_zeta[4 + gamma][alpha]*G_dual[alpha][beta]*scal4;
					}

					ierr = fem->dot_product(Ereal1,v1,&scal5); CHKERRQ(ierr);
					ierr = fem->dot_product(Ereal2,v1,&scal6); CHKERRQ(ierr);

					//term 4
					A11_res4[p] = A11_res4[p] + fem->dpsi[i][p][beta]*A_el_gp[gamma]*G_dual[alpha][beta]*scal5;
					A12_res4[p] = A12_res4[p] + fem->dpsi[i][p][beta]*A_el_gp[gamma]*G_dual[alpha][beta]*scal6;
					A21_res4[p] = A21_res4[p] + fem->dpsi[i][p][beta]*A_el_gp[2 + gamma]*G_dual[alpha][beta]*scal5;
					A22_res4[p] = A22_res4[p] + fem->dpsi[i][p][beta]*A_el_gp[2 + gamma]*G_dual[alpha][beta]*scal6;

					if (problem_type == forward){
						A31_res4[p] = A31_res4[p] + fem->dpsi[i][p][beta]*A_el_gp[4 + gamma]*G_dual[alpha][beta]*scal5;
						A32_res4[p] = A32_res4[p] + fem->dpsi[i][p][beta]*A_el_gp[4 + gamma]*G_dual[alpha][beta]*scal6;
					}

				}
			}
		}
	}

	return 0;
}

PetscErrorCode Dual_Solve::regularisation_stiffness_terms_calculation(int i, MyMat <double> E_real_el, MyMat <double> G_dual,
																	MyMat <double> dE1_zeta, MyMat <double> dE2_zeta,
																	MyVec <double> &A11_stiff1, MyVec <double> &A12_stiff1,
																	MyVec <double> &A21_stiff1, MyVec <double> &A22_stiff1,
																	MyVec <double> &A31_stiff1, MyVec <double> &A32_stiff1,
																	MyMat <double> &A11_stiff2, MyMat <double> &A12_stiff2,
																	MyMat <double> &A21_stiff2, MyMat <double> &A22_stiff2,
																	MyMat <double> &A31_stiff2, MyMat <double> &A32_stiff2,
																	MyMat <double> &A11_stiff3, MyMat <double> &A12_stiff3,
																	MyMat <double> &A21_stiff3, MyMat <double> &A22_stiff3,
																	MyMat <double> &A31_stiff3, MyMat <double> &A32_stiff3)
{

	PetscErrorCode ierr;

	MyVec <double> v1(grid->X_DOF);
	MyVec <double> v2(grid->X_DOF);
	MyVec <double> v3(grid->X_DOF);
	MyVec <double> v4(grid->X_DOF);

	MyVec <double> Ereal1(grid->X_DOF);
	MyVec <double> Ereal2(grid->X_DOF);

	double scal1, scal2, scal3, scal4;
	double scal5, scal6, scal7, scal8;
	double scal9, scal10, scal11, scal12;

	// E real
	Ereal1[0] = E_real_el[i][0];
	Ereal1[1] = E_real_el[i][1];
	Ereal1[2] = E_real_el[i][2];

	Ereal2[0] = E_real_el[i][grid->X_DOF];
	Ereal2[1] = E_real_el[i][grid->X_DOF + 1];
	Ereal2[2] = E_real_el[i][grid->X_DOF + 2];

	for (int k=0; k<grid->A_DOF; k++){

		A11_stiff1[k] = 0.0; A12_stiff1[k] = 0.0; A21_stiff1[k] = 0.0; 
		A22_stiff1[k] = 0.0; A31_stiff1[k] = 0.0; A32_stiff1[k] = 0.0;

		for (int p=0; p<fem->node_per_el; p++){

			A11_stiff2[k][p] = 0.0; A12_stiff2[k][p] = 0.0; A21_stiff2[k][p] = 0.0; 
			A22_stiff2[k][p] = 0.0; A31_stiff2[k][p] = 0.0; A32_stiff2[k][p] = 0.0;
			
			A11_stiff3[k][p] = 0.0; A12_stiff3[k][p] = 0.0; A21_stiff3[k][p] = 0.0; 
			A22_stiff3[k][p] = 0.0; A31_stiff3[k][p] = 0.0; A32_stiff3[k][p] = 0.0;

		}
	}

	for (int alpha=0; alpha<2; alpha++)
	{
		for (int beta=0; beta<2; beta++)
		{

			v1[0] = dE1_zeta[0][alpha];
			v1[1] = dE1_zeta[1][alpha];
			v1[2] = dE1_zeta[2][alpha];

			v2[0] = dE2_zeta[0][alpha];
			v2[1] = dE2_zeta[1][alpha];
			v2[2] = dE2_zeta[2][alpha];

			v3[0] = dE1_zeta[0][beta];
			v3[1] = dE1_zeta[1][beta];
			v3[2] = dE1_zeta[2][beta];

			v4[0] = dE2_zeta[0][beta];
			v4[1] = dE2_zeta[1][beta];
			v4[2] = dE2_zeta[2][beta];

			scal1 = 0.0; scal2 = 0.0; scal3 = 0.0; scal4 = 0.0; 
			scal5 = 0.0; scal6 = 0.0; scal7 = 0.0; scal8 = 0.0;
			scal9 = 0.0; scal10 = 0.0; scal11 = 0.0; scal12 = 0.0;
			
			ierr = fem->dot_product(v1,v3,&scal1); CHKERRQ(ierr);
			ierr = fem->dot_product(v2,v3,&scal2); CHKERRQ(ierr);
			
			ierr = fem->dot_product(v1,v4,&scal3); CHKERRQ(ierr);
			ierr = fem->dot_product(v2,v4,&scal4); CHKERRQ(ierr);

			ierr = fem->dot_product(Ereal1,v3,&scal5); CHKERRQ(ierr);
			ierr = fem->dot_product(Ereal2,v3,&scal6); CHKERRQ(ierr);

			ierr = fem->dot_product(Ereal1,v4,&scal7); CHKERRQ(ierr);
			ierr = fem->dot_product(Ereal2,v4,&scal8); CHKERRQ(ierr);

			ierr = fem->dot_product(Ereal1,v1,&scal9); CHKERRQ(ierr);
			ierr = fem->dot_product(Ereal1,v2,&scal10); CHKERRQ(ierr);

			ierr = fem->dot_product(Ereal2,v1,&scal11); CHKERRQ(ierr);
			ierr = fem->dot_product(Ereal2,v2,&scal12); CHKERRQ(ierr);

			// term 2
			A11_stiff1[0] = A11_stiff1[0] + scal1*G_dual[alpha][beta];
			A11_stiff1[1] = A11_stiff1[1] + scal2*G_dual[alpha][beta];

			A12_stiff1[0] = A12_stiff1[0] + scal3*G_dual[alpha][beta];
			A12_stiff1[1] = A12_stiff1[1] + scal4*G_dual[alpha][beta];

			A21_stiff1[2] = A21_stiff1[2] + scal1*G_dual[alpha][beta];
			A21_stiff1[3] = A21_stiff1[3] + scal2*G_dual[alpha][beta];

			A22_stiff1[2] = A22_stiff1[2] + scal3*G_dual[alpha][beta];
			A22_stiff1[3] = A22_stiff1[3] + scal4*G_dual[alpha][beta];

			if (problem_type == forward)
			{
				A31_stiff1[4] = A31_stiff1[4] + scal1*G_dual[alpha][beta];
				A31_stiff1[5] = A31_stiff1[5] + scal2*G_dual[alpha][beta];

				A32_stiff1[4] = A32_stiff1[4] + scal3*G_dual[alpha][beta];
				A32_stiff1[5] = A32_stiff1[5] + scal4*G_dual[alpha][beta];
			}

			for (int p=0;p<fem->node_per_el;p++)
			{

				// term 3
				A11_stiff2[0][p] = A11_stiff2[0][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal5;
				A11_stiff2[1][p] = A11_stiff2[1][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal6;

				A12_stiff2[0][p] = A12_stiff2[0][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal7;
				A12_stiff2[1][p] = A12_stiff2[1][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal8;

				A21_stiff2[2][p] = A21_stiff2[2][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal5;
				A21_stiff2[3][p] = A21_stiff2[3][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal6;

				A22_stiff2[2][p] = A22_stiff2[2][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal7;
				A22_stiff2[3][p] = A22_stiff2[3][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal8;

				if (problem_type == forward){

					A31_stiff2[4][p] = A31_stiff2[4][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal5;
					A31_stiff2[5][p] = A31_stiff2[5][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal6;

					A32_stiff2[4][p] = A32_stiff2[4][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal7;
					A32_stiff2[5][p] = A32_stiff2[5][p] + fem->dpsi[i][p][alpha]*G_dual[alpha][beta]*scal8;
				}

				// term 4
				A11_stiff3[0][p] = A11_stiff3[0][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal9;
				A11_stiff3[1][p] = A11_stiff3[1][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal10;

				A12_stiff3[0][p] = A12_stiff3[0][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal11;
				A12_stiff3[1][p] = A12_stiff3[1][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal12;

				A21_stiff3[2][p] = A21_stiff3[2][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal9;
				A21_stiff3[3][p] = A21_stiff3[3][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal10;

				A22_stiff3[2][p] = A22_stiff3[2][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal11;
				A22_stiff3[3][p] = A22_stiff3[3][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal12;

				if (problem_type == forward){

					A31_stiff3[4][p] = A31_stiff3[4][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal9;
					A31_stiff3[5][p] = A31_stiff3[5][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal10;

					A32_stiff3[4][p] = A32_stiff3[4][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal11;
					A32_stiff3[5][p] = A32_stiff3[5][p] + fem->dpsi[i][p][beta]*G_dual[alpha][beta]*scal12;
				}

			}
		}
	}
	
	
	return 0;
}

PetscErrorCode Dual_Solve::global_newton_solve_bcs(Mat &K_dual, Vec &rhs_vec, Vec &delta_dual_vec)
{

	PetscErrorCode ierr;

	int tbel = grid->nbel*grid->A_DOF;
  	int rows[grid->A_DOF];
  	double val[grid->A_DOF];

  	int rows_full[tbel];

	if (problem_type == inverse)
	{
		for (int i=0;i<grid->nbel;i++){

			// rows[0] = grid->dof_per_node*grid->b_el_conn[i][0];
			// rows[1] = grid->dof_per_node*grid->b_el_conn[i][0]+1;
			rows[0] = grid->dof_per_node*grid->b_el_conn[i][0]+2;
			rows[1] = grid->dof_per_node*grid->b_el_conn[i][0]+3;
			rows[2] = grid->dof_per_node*grid->b_el_conn[i][0]+4;
			rows[3] = grid->dof_per_node*grid->b_el_conn[i][0]+5;

			// val[0] = 0.0;
			// val[1] = 0.0;
			val[0] = 0.0;
			val[1] = 0.0;
			val[2] = 0.0;
			val[3] = 0.0;

			// rows_full[2*i] = rows[0];
			// rows_full[2*i+1] = rows[1];
			rows_full[grid->A_DOF*i] = rows[0];
			rows_full[grid->A_DOF*i+1] = rows[1];
			rows_full[grid->A_DOF*i+2] = rows[2];
			rows_full[grid->A_DOF*i+3] = rows[3];

			ierr = VecSetValues(delta_dual_vec,grid->A_DOF,rows,val,INSERT_VALUES); CHKERRQ(ierr);
			ierr = VecSetValues(rhs_vec,grid->A_DOF,rows,val,INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	else if (problem_type == forward)
	{
		for (int i=0;i<grid->nbel;i++){

			// rows[0] = grid->dof_per_node*grid->b_el_conn[i][0];
			// rows[1] = grid->dof_per_node*grid->b_el_conn[i][0]+1;
			rows[0] = grid->dof_per_node*grid->b_el_conn[i][0]+2;
			rows[1] = grid->dof_per_node*grid->b_el_conn[i][0]+3;
			rows[2] = grid->dof_per_node*grid->b_el_conn[i][0]+4;
			rows[3] = grid->dof_per_node*grid->b_el_conn[i][0]+5;
			rows[4] = grid->dof_per_node*grid->b_el_conn[i][0]+6;
			rows[5] = grid->dof_per_node*grid->b_el_conn[i][0]+7;

			// val[0] = 0.0;
			// val[1] = 0.0;
			val[0] = 0.0;
			val[1] = 0.0;
			val[2] = 0.0;
			val[3] = 0.0;
			val[4] = 0.0;
			val[5] = 0.0;

			// rows_full[2*i] = rows[0];
			// rows_full[2*i+1] = rows[1];
			rows_full[grid->A_DOF*i] = rows[0];
			rows_full[grid->A_DOF*i+1] = rows[1];
			rows_full[grid->A_DOF*i+2] = rows[2];
			rows_full[grid->A_DOF*i+3] = rows[3];
			rows_full[grid->A_DOF*i+4] = rows[4];
			rows_full[grid->A_DOF*i+5] = rows[5];

			ierr = VecSetValues(delta_dual_vec,grid->A_DOF,rows,val,INSERT_VALUES); CHKERRQ(ierr);
			ierr = VecSetValues(rhs_vec,grid->A_DOF,rows,val,INSERT_VALUES); CHKERRQ(ierr);
		}
	}

	// final assembly of global vector
	ierr = VecAssemblyBegin(rhs_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs_vec); CHKERRQ(ierr);

	// final assembly of global vector
	ierr = VecAssemblyBegin(delta_dual_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(delta_dual_vec); CHKERRQ(ierr);
	
	ierr = MatZeroRows(K_dual,tbel,rows_full,1.0,delta_dual_vec,rhs_vec); CHKERRQ(ierr);

	return 0;

}


PetscErrorCode Dual_Solve::get_dual_l2_norm()
{

	PetscErrorCode ierr;

	dual_S_total = 0.0; // initialised to zero

	qf_global = 0.0; // initialised to zero

	// undef_global = 0.0; // initialised to zero

	double S_m1_el = 0.0 , S_m2_el = 0.0, S_F_dx_el = 0.0, S_quad_H_el = 0.0;

	double qfactor_el = 0.0, compatibility_el = 0.0;

	double qf_g = 0.0, undef_g = 0.0;

	MyMat<double> beta_el(fem->node_per_el,MyVec<double> (grid->BETA_DOF));

	MyMat<double> A_el(fem->node_per_el, MyVec<double>(grid->A_DOF));

	MyMat <double> x_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> x_0_gp_el(fem->ngp,MyVec<double>(this->dim));

	MyMat <double> F_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyMat <double> F_0_dual_gp_el(fem->ngp,MyVec<double>(F_DOF));

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	for (int ie=0;ie<grid->nel;++ie)
	{	

		// getting dual fields for each elem  
		for (int i=0;i<fem->node_per_el;i++){
			
			for (int k=0;k<grid->BETA_DOF;k++){
				beta_el[i][k] = beta[grid->el_conn[ie][i]][k];
			}

			for (int k=0;k<grid->A_DOF;k++){
				A_el[i][k] = A_dual[grid->el_conn[ie][i]][k];
			}
		}

		// getting basis vectors and jacobian for each elem
		for (int i=0;i<fem->ngp;i++){
			for (int k=0;k<2*grid->X_DOF;k++){
				E_dual_el[i][k] = fem->E_dual[ie][i][k];
				//std::cout << "E_dual_el for i and k is " << E_dual_el[i][k] << std::endl;
			}
			
			detJ_e[i] = fem->detJ_G[ie][i];	

		}

		// storing primal fields at gpts
		for (int i=0;i<fem->ngp;i++){

			for (int k=0;k<dim;k++){
				x_gp_el[i][k] = x_gp[ie][i][k];
			}

			for (int k=0;k<dim;k++){
				x_0_gp_el[i][k] = x_0_gp[ie][i][k];
			}

			for (int k=0;k<F_DOF;k++){
				F_gp_el[i][k] = F_gp[ie][i][k];
			}

			for (int k=0;k<F_DOF;k++){
				F_0_dual_gp_el[i][k] = F_0_dual_gp[ie][i][k];
			}
		}

		// elem dual l2 norm
		ierr = this->get_elem_dual_l2_norm(beta_el,A_el,E_dual_el,detJ_e,x_gp_el,x_0_gp_el,F_gp_el,F_0_dual_gp_el,
												&S_m1_el, &S_m2_el, &S_F_dx_el, &S_quad_H_el,
												&qfactor_el, &compatibility_el,&qf_g,&undef_g); CHKERRQ(ierr);

		dual_S_mu1[ie] = S_m1_el;
		dual_S_mu2[ie] = S_m2_el;
		dual_S_F_dx[ie] = S_F_dx_el;
		dual_S_quad_H[ie] = S_quad_H_el;
		qfactor[ie] = qfactor_el;
		compatibility_check[ie] = compatibility_el;

		dual_S_norm[ie] = S_m1_el + S_m2_el + S_F_dx_el + S_quad_H_el;

		dual_S_total = dual_S_total +  dual_S_norm[ie];

		qf_global = qf_global + qf_g;
		// undef_global = undef_global + undef_g;

	} 

	return 0;

}


PetscErrorCode Dual_Solve::get_elem_dual_l2_norm(MyMat <double> beta_el, MyMat <double> A_el, 
												 MyMat <double> E_dual_el,  MyVec<double> detJ_e, 
												 MyMat <double> x_gp_el, MyMat <double> x_0_gp_el,
												 MyMat <double> F_gp_el, MyMat <double> F_0_dual_gp_el, 
												 double *S_m1_el, double *S_m2_el, 
												 double *S_F_dx_el, double *S_quad_H_el,
												 double *qfactor_el, double *compatibility_el,
												 double *qf_g, double *undef_g)
{

	PetscErrorCode ierr;

	MyMat<double> F(dim,MyVec<double>(2));

	MyMat<double> F_dual(dim,MyVec<double>(2));

	MyMat<double> F_0(dim,MyVec<double>(2));

	MyMat<double> F_0_dual(dim,MyVec<double>(2));

	MyMat<double> G_dual(2,MyVec<double>(2));

	MyMat<double> G_real(2,MyVec<double>(2));

	MyMat<double> C(2,MyVec<double>(2));

	double tr_C, det_C;

	double mu1, mu2;

	double det_F;

	int index;

	double S_m1, S_m2, S_F_dx, S_quad_H;

	double beta_1_gp, beta_2_gp, A_at_gp;

	double undef_area = 0.0;

	double qf, compatibility;

	*S_m1_el = 0.0;
	*S_m2_el = 0.0;
	*S_F_dx_el = 0.0;
	*S_quad_H_el = 0.0;
	*qfactor_el = 0.0;
	*compatibility_el = 0.0;
	*qf_g = 0.0;
	*undef_g = 0.0;

	MyMat <double> dx_zeta(dim,MyVec<double>(2));

	// loop over the gauss pts
	for (int i=0;i<fem->ngp;i++)
	{	

		S_F_dx = 0.0;
		S_quad_H = 0.0;
		compatibility = 0.0;

		for (int k1=0;k1<dim;k1++){
			for(int k2=0;k2<2;k2++){

				index = k1*2 + k2;
				F[k1][k2] = F_gp_el[i][index];
				F_0_dual[k1][k2] = F_0_dual_gp_el[i][index];
			}
		}

		ierr = this->get_G_dual(i,E_dual_el,G_dual); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);		

		ierr = this->tranform_F_to_dual(F,G_dual,F_dual); CHKERRQ(ierr);

		ierr = this->tranform_F_from_dual(F_0_dual,G_real,F_0); CHKERRQ(ierr);

		ierr = this->get_C(F,F_dual,C); CHKERRQ(ierr);

		ierr = this->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = this->determinant_of_C(F,G_real,&det_C); CHKERRQ(ierr);
		
		ierr = this->get_invariants(tr_C,det_C,&mu1,&mu2); CHKERRQ(ierr);

		ierr = this->get_beta_field_at_gp(beta_el,i,&beta_1_gp,&beta_2_gp); CHKERRQ(ierr);


		// S_m1 = std::pow(std::abs(mu1-lambda_1*lambda_1),1.0+EPS)*beta_1_gp;
		// S_m2 = std::pow(std::abs(mu2-lambda_2*lambda_2),1.0+EPS)*beta_2_gp;
		S_m1 = (mu1-(lambda_1*lambda_1))*beta_1_gp;
		S_m2 = (mu2-(lambda_2*lambda_2))*beta_2_gp;
		qf = (mu1- (lambda_1*lambda_1))*(mu1- (lambda_1*lambda_1)) + 
				(mu2- (lambda_2*lambda_2))*(mu2- (lambda_2*lambda_2));

		for (int d=0;d<this->dim;d++){
		
			for (int alpha=0;alpha<2;alpha++){

				dx_zeta[d][alpha] = 0.0; // initializing to zero

				for (int k=0;k<fem->ngp;k++){

					dx_zeta[d][alpha] = dx_zeta[d][alpha] + fem->dpsi[i][k][alpha]*x_gp_el[k][d];
				}

			}
		}

		// check for det F (throws error if det (dx_zeta) < 1e-4)
		ierr = this->determinant_of_F(dx_zeta,G_real,&det_F); CHKERRQ(ierr);


		for (int k1=0;k1<this->dim; k1++){

			for (int k2=0;k2<2;k2++){

				index = k1*2 + k2;

				ierr = this->get_A_field_at_gp(A_el,index,i,&A_at_gp); CHKERRQ(ierr);

				S_F_dx = S_F_dx + A_at_gp*(F[k1][k2] - dx_zeta[k1][k2]);

				compatibility = compatibility + (F[k1][k2] - dx_zeta[k1][k2])*(F[k1][k2] - dx_zeta[k1][k2]);

				S_quad_H = S_quad_H + 0.5*this->pf*(F[k1][k2]-F_0[k1][k2])*(F_dual[k1][k2]-F_0_dual[k1][k2]);
					
			}
		}

		if (problem_type == inverse)
		{
			S_quad_H  = S_quad_H + 0.5*this->px*((x_gp_el[i][0] - x_0_gp_el[i][0])*(x_gp_el[i][0] - x_0_gp_el[i][0]) 
											+ (x_gp_el[i][1] - x_0_gp_el[i][1])*(x_gp_el[i][1] - x_0_gp_el[i][1]));
		}
		else if (problem_type == forward)
		{
			S_quad_H  = S_quad_H + 0.5*this->px*((x_gp_el[i][0] - x_0_gp_el[i][0])*(x_gp_el[i][0] - x_0_gp_el[i][0]) 
											+ (x_gp_el[i][1] - x_0_gp_el[i][1])*(x_gp_el[i][1] - x_0_gp_el[i][1])
											+ (x_gp_el[i][2] - x_0_gp_el[i][2])*(x_gp_el[i][2] - x_0_gp_el[i][2]));
		}


		*S_m1_el = *S_m1_el + S_m1*detJ_e[i]*fem->wt[i];
		*S_m2_el = *S_m2_el + S_m2*detJ_e[i]*fem->wt[i];
		*S_F_dx_el = *S_F_dx_el + S_F_dx*detJ_e[i]*fem->wt[i];
		*S_quad_H_el = *S_quad_H_el + S_quad_H*detJ_e[i]*fem->wt[i];

		*qfactor_el = *qfactor_el + qf*detJ_e[i]*fem->wt[i];
		*compatibility_el = *compatibility_el + compatibility*detJ_e[i]*fem->wt[i];

		undef_area = undef_area + detJ_e[i]*fem->wt[i];

	}// loop over gpts

	*S_m1_el = *S_m1_el/undef_area;
	*S_m2_el = *S_m2_el/undef_area;
	*S_F_dx_el = *S_F_dx_el/undef_area;
	*S_quad_H_el = *S_quad_H_el/undef_area;

	*qf_g = *qfactor_el;
	*undef_g = undef_area;

	*qfactor_el = std::sqrt(*qfactor_el/undef_area);
	*compatibility_el = std::sqrt(*compatibility_el/undef_area);

	return 0;
}

PetscErrorCode Dual_Solve::get_derivative_of_mus(MyMat <double> C_inv, double tr_C, double det_C,
												MyMat <double> G_dual, MyMat <double> F_guess, 
												int k1, int k2, double *d_mu1_F, double *d_mu2_F)
{

	double tmp = tr_C*tr_C - 4.0*det_C;

	double sqrt_tmp = 0.0;

	*d_mu1_F = 0.0;
	*d_mu2_F = 0.0;

	double term11 = 0.0, term12 = 0.0;
	double term21 = 0.0, term22 = 0.0;

	if (tmp > 1e-10){
		sqrt_tmp = std::sqrt(tmp);
	}
	else{
		sqrt_tmp = 0.0;
		throw "eigenvalues are same";
	}

	double one_by_sqrt_tmp;

	if (sqrt_tmp > 1e-5){
		one_by_sqrt_tmp = 1.0/sqrt_tmp;
	}
	else{
		one_by_sqrt_tmp = 0.0;
		// std::cout << "sqrt of tr_C^2 - 4*det_C is " << sqrt_tmp << std::endl;
		throw "divide by zero in derivatives of mu's calculation";
		// one_by_sqrt_tmp = 0.0;
	}


	for (int alpha=0;alpha<2;alpha++)
	{
		for (int beta=0;beta<2;beta++)
		{
			for (int m=0;m<this->dim;m++)
			{

				term11 = term11 + 0.5*(Id[beta][alpha] + 0.5*one_by_sqrt_tmp*(2.0*tr_C*Id[beta][alpha]
									- 4.0*det_C*C_inv[beta][alpha]))*(II[k1][m]*G_dual[k2][beta]*F_guess[m][alpha]);

				term21 = term21 + 0.5*(Id[beta][alpha] - 0.5*one_by_sqrt_tmp*(2.0*tr_C*Id[beta][alpha]
								- 4.0*det_C*C_inv[beta][alpha]))*(II[k1][m]*G_dual[k2][beta]*F_guess[m][alpha]);

				for (int rho=0;rho<2;rho++)
				{

					term12 = term12 + 0.5*(Id[beta][alpha] + 0.5*one_by_sqrt_tmp*(2.0*tr_C*Id[beta][alpha]
									- 4.0*det_C*C_inv[beta][alpha]))*(II[m][k1]*G_dual[rho][beta]*F_guess[m][rho]*Id[alpha][k2]);

					term22 = term22 + 0.5*(Id[beta][alpha] - 0.5*one_by_sqrt_tmp*(2.0*tr_C*Id[beta][alpha]
									- 4.0*det_C*C_inv[beta][alpha]))*(II[m][k1]*G_dual[rho][beta]*F_guess[m][rho]*Id[alpha][k2]);
				}
			}
		}
	}

	*d_mu1_F = term11 + term12;
	*d_mu2_F = term21 + term22;

	return 0;
}


PetscErrorCode Dual_Solve::first_derivative_of_invariants(MyVec <double> eig_v1, MyVec <double> eig_v2,
															MyMat <double> G_dual, MyMat<double> F_guess,
															MyVec <double> &d_mu1_F, MyVec <double> &d_mu2_F)
{	

	int index;

	// zeroing the output first as summation in formula
	for (int k1=0;k1<F_DOF; k1++){

		d_mu1_F[k1] = 0.0;
		d_mu2_F[k1] = 0.0;
	}

	MyMat <double> d_mu1_C (2,MyVec <double> (2));
	MyMat <double> d_mu2_C (2,MyVec <double> (2));

	MyVec <double> term11(F_DOF);
	MyVec <double> term12(F_DOF);
	MyVec <double> term21(F_DOF);
	MyVec <double> term22(F_DOF);


	for (int alpha=0;alpha<2;alpha++){

		for (int beta=0;beta<2;beta++){

			for (int rho=0;rho<2;rho++){

				d_mu1_C[alpha][beta] =  d_mu1_C[alpha][beta] + eig_v1[rho]*G_dual[rho][alpha]*eig_v1[beta];
				d_mu2_C[alpha][beta] =  d_mu2_C[alpha][beta] + eig_v2[rho]*G_dual[rho][alpha]*eig_v2[beta];
			}
			
		}
	}


	for (int k1=0;k1<this->dim; k1++)
	{	
		for (int k2=0; k2<2; k2++)
		{	

			index = k1*2 + k2;

			for (int alpha=0;alpha<2;alpha++)
			{
				for (int beta=0;beta<2;beta++)
				{
					for (int m=0;m<this->dim;m++)
					{	

						term11[index] = term11[index] + d_mu1_C[alpha][beta]*(II[k1][m]*G_dual[k2][beta]*F_guess[m][alpha]);

						term21[index] = term21[index] + d_mu2_C[alpha][beta]*(II[k1][m]*G_dual[k2][beta]*F_guess[m][alpha]);

						for (int rho=0;rho<2;rho++)
						{

							term12[index] = term12[index] + d_mu1_C[alpha][beta]*(II[m][k1]*G_dual[rho][beta]*F_guess[m][rho]*Id[alpha][k2]);

							term22[index] = term22[index] + d_mu2_C[alpha][beta]*(II[m][k1]*G_dual[rho][beta]*F_guess[m][rho]*Id[alpha][k2]);
						}
					}
				}
			}
		}
	}


	for (int k1=0;k1<this->dim; k1++)
	{	
		for (int k2=0; k2<2; k2++)
		{	

			index = k1*2 + k2;

			d_mu1_F[index] = d_mu1_F[index] + term11[index] + term12[index];
			d_mu2_F[index] = d_mu2_F[index] + term21[index] + term22[index];

		}
	}

	return 0;
}

PetscErrorCode Dual_Solve::second_derivative_of_invariants(MyVec <double> eig_v1,MyVec <double> eig_v2,
															double mu1, double mu2, double tr_C, MyMat <double> C,
															MyMat <double> G_dual, MyMat <double> G_real, MyMat <double> F_guess,
															MyMat <double> &dd_mu1_FF, MyMat <double> &dd_mu2_FF)
{	

	int index1, index2;

	int ind1, ind2;

	double tmp_pi_theta = 0.0;

	double scal_theta_p  = 0.0;
	double scal_pi_p = 0.0;
	double scal_mu_k = 0.0;
	double scal_rho_k = 0.0;

	// zeroing the output first as summation in formula
	for (int k1=0;k1<F_DOF; k1++){
		for (int k2=0; k2<F_DOF; k2++){

			dd_mu1_FF[k1][k2] = 0.0;
			dd_mu2_FF[k1][k2] = 0.0;	
		}
	}

	MyMat <double> dd_mu1_CC (4,MyVec <double> (4));
	MyMat <double> dd_mu2_CC (4,MyVec <double> (4));

	if (std::abs(2.0*mu1 - tr_C) > 1e-3){		

		for (int rho=0;rho<2;rho++){
			for (int mu=0;mu<2;mu++){

				index1 = 2*rho + mu;

				for (int pi=0;pi<2;pi++){
					for (int theta=0;theta<2;theta++){

						index2 = 2*pi + theta;

						tmp_pi_theta = 0.0;

						for (int alpha=0;alpha<2;alpha++){
							for (int beta=0;beta<2;beta++){

								tmp_pi_theta = tmp_pi_theta + C[alpha][beta]*G_dual[alpha][pi]*G_real[beta][theta];
							}
						}


						dd_mu1_CC[index1][index2] = dd_mu1_CC[index1][index2] + (1.0/(2.0*mu1 - tr_C))*(G_dual[rho][pi]*G_real[mu][theta]
																-(Id[rho][mu]*Id[pi][theta]
																- (eig_v1[0]*G_dual[0][rho] + eig_v1[1]*G_dual[1][rho])
																	*eig_v1[mu]*Id[pi][theta]))
																- (1.0/((2.0*mu1 - tr_C)*(2.0*mu1 - tr_C)))*(
																	(2.0*(eig_v1[0]*G_dual[0][rho] + eig_v1[1]*G_dual[1][rho])*
																	eig_v1[mu]) - Id[rho][mu])*(tmp_pi_theta - (tr_C - mu1)*Id[pi][theta]);

						dd_mu2_CC[index1][index2] = dd_mu2_CC[index1][index2] + (1.0/(2.0*mu2 - tr_C))*(G_dual[rho][pi]*G_real[mu][theta]
																-(Id[rho][mu]*Id[pi][theta]
																- (eig_v2[0]*G_dual[0][rho] + eig_v2[1]*G_dual[1][rho])
																	*eig_v2[mu]*Id[pi][theta]))
																- (1.0/((2.0*mu2 - tr_C)*(2.0*mu2 - tr_C)))*(
																	(2.0*(eig_v2[0]*G_dual[0][rho] + eig_v2[1]*G_dual[1][rho])*
																	eig_v2[mu]) - Id[rho][mu])*(tmp_pi_theta - (tr_C - mu2)*Id[pi][theta]); 											 

						
					}
				}
			}
		}
	}
	else 
	{
		// std::cout << "two invariants of C are same that is " << mu1 << " and " << mu2 << std::endl;
		// throw "divide by zero in second derivatives of mu's calculation";

		for (int rho=0;rho<2;rho++){
			for (int mu=0;mu<2;mu++){

				index1 = 2*rho + mu;

				for (int pi=0;pi<2;pi++){
					for (int theta=0;theta<2;theta++){

						index2 = 2*pi + theta;

						dd_mu1_CC[index1][index2] = dd_mu1_CC[index1][index2] + 0.5*(1.0/(std::sqrt(mu1)+  std::sqrt(mu2)))*
																				(eig_v1[0]*G_dual[0][rho] + 
																				eig_v1[1]*G_dual[1][rho])*eig_v2[mu]*(eig_v2[theta]*
																				(eig_v1[0]*G_dual[0][pi] + eig_v1[1]*G_dual[1][pi])
																				+ 	eig_v1[theta]*(eig_v2[0]*G_dual[0][pi]
																				 + eig_v2[1]*G_dual[1][pi]));

						dd_mu2_CC[index1][index2] = dd_mu2_CC[index1][index2] + 0.5*(1.0/(std::sqrt(mu1)+  std::sqrt(mu2)))*
																				(eig_v1[0]*G_dual[0][rho] + 
																				eig_v1[1]*G_dual[1][rho])*eig_v2[mu]*(eig_v2[theta]*
																				(eig_v1[0]*G_dual[0][pi] + eig_v1[1]*G_dual[1][pi])
																				+ 	eig_v1[theta]*(eig_v2[0]*G_dual[0][pi]
																				 + eig_v2[1]*G_dual[1][pi]));											 
						
					}
				}
			}
		}
	}


	for (int k=0;k<this->dim;k++){
		for (int gamma=0;gamma<2;gamma++){

			ind1 = k*2 + gamma;

			for (int p=0;p<this->dim;p++){
				for (int omega=0;omega<2;omega++){


					ind2 = p*2 + omega;


					for (int rho=0;rho<2;rho++){
						for (int mu=0;mu<2;mu++){

							index1 = rho*2 + mu;

							for (int pi=0;pi<2;pi++){
								for (int theta=0;theta<2;theta++){

									index2 = pi*2 + theta;	

									scal_theta_p  = 0.0;
									scal_pi_p = 0.0;
									scal_mu_k = 0.0;
									scal_rho_k = 0.0;

									for (int m=0; m<this->dim; m++){

										scal_pi_p = scal_pi_p + II[p][m]*F_guess[m][pi];
										scal_rho_k = scal_rho_k +  II[k][m]*F_guess[m][rho];

										for (int delta=0;delta<2;delta++){

											scal_theta_p = scal_theta_p +  II[m][p]*G_dual[delta][theta]*F_guess[m][delta];	
											scal_mu_k = scal_mu_k +  II[m][k]*G_dual[delta][mu]*F_guess[m][delta];
										}
									}

									dd_mu1_FF[ind1][ind2] = dd_mu1_FF[ind1][ind2] + dd_mu1_CC[index1][index2]*
																(scal_theta_p*Id[pi][omega] + scal_pi_p*G_dual[omega][theta])*
																(scal_mu_k*Id[rho][gamma] + scal_rho_k*G_dual[gamma][mu]);

									dd_mu2_FF[ind1][ind2] = dd_mu2_FF[ind1][ind2] + dd_mu2_CC[index1][index2]*
																(scal_theta_p*Id[pi][omega] + scal_pi_p*G_dual[omega][theta])*
																(scal_mu_k*Id[rho][gamma] + scal_rho_k*G_dual[gamma][mu]);																
								}
							}
						}
					}


				}
			}
		}
	}




	return 0;
}



inline PetscErrorCode Dual_Solve::tranform_C_on_both_legs(MyMat <double> G_real,MyMat <double> G_dual,
														MyMat <double> C, MyMat <double> &C_tranf)
{ 

	// C_tranf = C_tranf^{\alpha}_{\beta} E_{\alpha} \otimes E^{\beta}

	for (int alpha=0;alpha<2;alpha++){
		for (int beta=0;beta<2;beta++){

			C_tranf[alpha][beta] = 0.0;

			for (int gamma=0;gamma<2;gamma++){
				for (int delta=0;delta<2;delta++){

					C_tranf[alpha][beta] = C_tranf[alpha][beta] + 
												G_dual[gamma][alpha]*C[gamma][delta]*G_real[delta][beta];
				}
			}
		}
	}

	return 0;
}												


inline PetscErrorCode Dual_Solve::get_G_dual(int gp, MyMat <double> E_dual_el, MyMat<double> &G_dual)
{	

	PetscErrorCode ierr;

	MyVec <double> v1(3), v2(3);

	double val = 0.0;

	for (int k1=0;k1<2;k1++){
		for (int k2=0;k2<2;k2++){

			v1[0] = E_dual_el[gp][3*k1];
			v1[1] = E_dual_el[gp][3*k1+1];
			v1[2] = E_dual_el[gp][3*k1+2];

			v2[0] = E_dual_el[gp][3*k2];
			v2[1] = E_dual_el[gp][3*k2+1];
			v2[2] = E_dual_el[gp][3*k2+2];


			ierr = fem->dot_product(v1,v2,&val); CHKERRQ(ierr) ;
			G_dual[k1][k2] = val;

		}
	}

	return 0;

}


inline PetscErrorCode Dual_Solve::get_H_term(MyMat <double> F_guess, 
											MyMat <double> G_dual, int k1, int k2,
											double *h_term)
{
 	
 	*h_term = 0.0;

	for (int m=0;m<dim;m++){
		for (int delta=0;delta<2;delta++){

			*h_term = *h_term + F_guess[m][delta]*II[k1][m]*G_dual[k2][delta];
		}
	}

	// std::cout << "h term for k1 and k2 " << k1 << " and " << k2 << " is "  << *h_term << std::endl;


	return 0;
}											

inline PetscErrorCode Dual_Solve::get_A_field_at_gp(MyMat <double> A_el,
													int index, int gp, double *A_at_gp)
{	

	*A_at_gp = 0.0;

	for (int i=0;i<fem->node_per_el;i++)
	{
		*A_at_gp = *A_at_gp + fem->psi[gp][i]*A_el[i][index];

	}

	return 0;
}


inline PetscErrorCode Dual_Solve::get_beta_field_at_gp(MyMat <double> beta_el, int gp, 
														double *beta_1_gp, double *beta_2_gp)
{	

	*beta_1_gp = 0.0;
	*beta_2_gp = 0.0;

	for (int i=0;i<fem->node_per_el;i++){

		*beta_1_gp = *beta_1_gp + fem->psi[gp][i]*beta_el[i][0];
		*beta_2_gp = *beta_2_gp + fem->psi[gp][i]*beta_el[i][1];

	}

	return 0;
}

inline PetscErrorCode Dual_Solve::get_derivative_of_beta_at_gp(MyMat <double> beta_el,int gp,
															MyVec <double> beta_1_zeta,
															MyVec <double> beta_2_zeta)
{

	for (int k=0;k<2;k++){
		beta_1_zeta[k] = 0.0;
		beta_2_zeta[k] = 0.0;
	}

	for (int i=0;i<fem->node_per_el;i++){

		beta_1_zeta[0] = beta_1_zeta[0] + fem->dpsi[gp][i][0]*beta_el[i][0];
		beta_1_zeta[1] = beta_1_zeta[1] + fem->dpsi[gp][i][1]*beta_el[i][0];

		beta_2_zeta[0] = beta_2_zeta[0] + fem->dpsi[gp][i][0]*beta_el[i][1];
		beta_2_zeta[1] = beta_2_zeta[1] + fem->dpsi[gp][i][1]*beta_el[i][1];
	}

	return 0;
}

inline PetscErrorCode Dual_Solve::get_div_A_field_at_gp(MyMat<double> A_el, int gp, double *div_A)
{
	if (problem_type == forward){
		throw "get_div_A_field_at_gp function need to made dim independent";
	}

	double d_A11_zeta1, d_A12_zeta2, d_A21_zeta1, d_A22_zeta2;

	d_A11_zeta1 = A_el[0][0]*fem->dpsi[gp][0][0] + A_el[1][0]*fem->dpsi[gp][1][0] + 
						A_el[2][0]*fem->dpsi[gp][2][0] +  A_el[3][0]*fem->dpsi[gp][3][0];
	
	d_A12_zeta2 = A_el[0][1]*fem->dpsi[gp][0][1] + A_el[1][1]*fem->dpsi[gp][1][1] + 
						A_el[2][1]*fem->dpsi[gp][2][1] +  A_el[3][1]*fem->dpsi[gp][3][1];

	d_A21_zeta1 = A_el[0][2]*fem->dpsi[gp][0][0] + A_el[1][2]*fem->dpsi[gp][1][0] + 
						A_el[2][2]*fem->dpsi[gp][2][0] +  A_el[3][2]*fem->dpsi[gp][3][0];
	
	d_A22_zeta2 = A_el[0][3]*fem->dpsi[gp][0][1] + A_el[1][3]*fem->dpsi[gp][1][1] + 
						A_el[2][3]*fem->dpsi[gp][2][1] +  A_el[3][3]*fem->dpsi[gp][3][1];

	div_A[0] = d_A11_zeta1 + d_A12_zeta2;
	div_A[1] = d_A21_zeta1 + d_A22_zeta2;

	return 0;
}

inline PetscErrorCode Dual_Solve::get_invariants(double tr_C, double det_C, double *mu1, double *mu2)
{	

	double tmp = std::pow(tr_C,2) - 4.0*det_C;

	double tmp2;

	if (tmp>0.0)
		tmp2 = std::sqrt(tmp);
	else
		tmp2 = 0.0;
	

	*mu1 = 0.5*(tr_C + tmp2);
	*mu2 = 0.5*(tr_C - tmp2);

	return 0;
}


inline PetscErrorCode Dual_Solve::get_eigenvectors_of_C(MyMat <double> C, double mu1, double mu2, MyMat <double> G_dual,
														MyVec <double> &eig_v1, MyVec <double> &eig_v2)
{

	double tmp, tmp2; 

	if (std::abs(C[0][0] - mu1) > 1e-4){

		eig_v1[0] = - C[0][1]/(C[0][0] - mu1);
		eig_v1[1] = 1.0;
	}
	else {
		
		if (std::abs(C[1][1] - mu1) > 1e-4){

			eig_v1[1] = - C[1][0]/(C[1][1] - mu1);
			eig_v1[0] = 1.0;

		}
		else {

			// std::cout << " division by zero in eigenvector calculation " << C[0][0] << " and " << mu1 << std::endl;
			// std::cout << "mu1 and mu2 are " << mu1 << " and " << mu2 << std::endl;
			// for (int i=0;i<2;i++){
			// 	for (int j=0;j<2;j++){
			// 		std::cout << "C[i][j] for i " << i << " and j " << j << " is " << C[i][j] << std::endl;
			// 	}
			// }
			// throw "division by zero in eigenvector calculation";
			eig_v1[1] = 0.0;
			eig_v1[0] = 1.0;

		}
	}

	// eig_v2[1] = 1.0;

	// if (std::abs(C[0][0] - mu2) > 1e-4){
	// 	eig_v2[0] = - C[0][1]/(C[0][0] - mu2);
	// }
	// else{
	// 	std::cout << " division by zero in eigenvector calculation " << C[0][0] << " and " << mu2 << std::endl;
	// 	std::cout << "mu1 and mu2 are " << mu1 << " and " << mu2 << std::endl;
	// 	for (int i=0;i<2;i++){
	// 		for (int j=0;j<2;j++){
	// 			std::cout << "C[i][j] for i " << i << " and j " << j << " is " << C[i][j] << std::endl;
	// 		}
	// 	}
	// 	throw "division by zero in eigenvector calculation";
	// }

	tmp = std::sqrt(eig_v1[0]*eig_v1[0]*G_dual[0][0] + eig_v1[1]*eig_v1[1]*G_dual[1][1] + 2.0*eig_v1[0]*eig_v1[1]*G_dual[0][1]);

	// tmp2 = std::sqrt(eig_v2[0]*eig_v2[0]*G_dual[0][0] + eig_v2[1]*eig_v2[1]*G_dual[1][1] + 2.0*eig_v2[0]*eig_v2[1]*G_dual[0][1]);

	for (int alpha=0;alpha<2;alpha++){

		eig_v1[alpha] = eig_v1[alpha]/tmp;
		// eig_v2[alpha] = eig_v2[alpha]/tmp2;

	}

	if (std::abs(eig_v1[1]*G_dual[1][1] + eig_v1[0]*G_dual[0][1]) < 1e-3){

		tmp2 = -(eig_v1[1]*G_dual[1][1] + eig_v1[0]*G_dual[0][1])/(eig_v1[0]*G_dual[0][0] + eig_v1[1]*G_dual[0][1]);

		eig_v2[1] = 1.0/(std::sqrt(tmp2*tmp2*G_dual[0][0] + G_dual[1][1] + 2.0*tmp2*G_dual[0][1]));
		eig_v2[0] = tmp2*eig_v2[1];
	}
	else {

		tmp2 = -(eig_v1[0]*G_dual[0][0] + eig_v1[1]*G_dual[0][1])/(eig_v1[1]*G_dual[1][1] + eig_v1[0]*G_dual[0][1]);

		eig_v2[0] = 1.0/(std::sqrt(G_dual[0][0] + tmp2*tmp2*G_dual[1][1] + 2.0*tmp2*G_dual[0][1]));
		eig_v2[1] = tmp2*eig_v2[0];

	}

	// std::cout << "eig_v1 is " << eig_v1[0] << " and " << eig_v1[1] << std::endl;
	// std::cout << "eig_v2 is " << eig_v2[0] << " and " << eig_v2[1] << std::endl;
	// std::cout << "tmp2 is " << tmp2 << std::endl;
	// std::cout << "scal1 is " << (eig_v1[0]*G_dual[0][0] + eig_v1[1]*G_dual[0][1]) << std::endl;
	// std::cout << "scal2 is " << (eig_v1[1]*G_dual[1][1] + eig_v1[0]*G_dual[0][1]) << std::endl;


	return 0;

}

inline PetscErrorCode Dual_Solve::trace_of_C(MyMat <double> C, double* tr_C)
{
	*tr_C = C[0][0] + C[1][1];

	return 0;
}

inline PetscErrorCode Dual_Solve::determinant_of_C(MyMat <double> F, MyMat <double> G, double *det_C)
{	

	PetscErrorCode ierr;
	
	double square_of_det_f, det_G;

	MyVec<double> v1(this->dim);
	MyVec<double> v2(this->dim);
	MyVec<double> v3(3);

	for (int i=0;i<this->dim;i++){
		v1[i] = F[i][0];
		v2[i] = F[i][1];		 
	}

	ierr = fem->cross_product(v1,v2,v3); CHKERRQ(ierr);

	// | \partial_1 x \cross \partial_2 x |^2
	square_of_det_f = v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];

	det_G = G[0][0]*G[1][1] - G[0][1]*G[1][0];
	
	*det_C = square_of_det_f/det_G;

	if (*det_C < 1e-10){
		std::cout << "determinant of the matrix C is " << *det_C << std::endl;
		throw "determinant of matrix is 0, not invertible";
	}

	return 0;
}


inline PetscErrorCode Dual_Solve::determinant_of_F(MyMat <double> dx_zeta, MyMat <double> G, double *det_F)
{	

	PetscErrorCode ierr;
	
	double det_f, sqrt_of_det_G;

	MyVec<double> v1(this->dim);
	MyVec<double> v2(this->dim);
	MyVec<double> v3(3);

	for (int i=0;i<this->dim;i++){
		v1[i] = dx_zeta[i][0];
		v2[i] = dx_zeta[i][1];		 
	}

	ierr = fem->cross_product(v1,v2,v3); CHKERRQ(ierr);

	// | \partial_1 x \cross \partial_2 x |
	det_f = std::sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);

	sqrt_of_det_G = std::sqrt(G[0][0]*G[1][1] - G[0][1]*G[1][0]);
	
	*det_F = det_f/sqrt_of_det_G;

	if (*det_F < 1e-6){
		std::cout << "determinant of the matrix dx_zeta is " << *det_F << std::endl;
		throw "determinant of matrix is 0, not invertible";
	}

	return 0;
}

PetscErrorCode Dual_Solve::tranform_F_to_dual(MyMat<double> F_gp, MyMat<double> G_dual, MyMat<double> &F_gp_dual)
{	

	// F_gp -->  F^k_\gamma c_k \otimes E^\gamma
	// F_gp_dual --> F_k^\gamma c^k \otimes E_\gamma

	for (int k=0;k<this->dim;k++){
		for (int gamma=0;gamma<2;gamma++){

			F_gp_dual[k][gamma] = 0.0;

			for (int p=0;p<this->dim;p++){
				for (int beta=0;beta<2;beta++){

					F_gp_dual[k][gamma] = F_gp_dual[k][gamma] + II[p][k]*G_dual[beta][gamma]*F_gp[p][beta];
				}
			}
		}
	}

	return 0;

}

PetscErrorCode Dual_Solve::tranform_F_from_dual(MyMat<double> F_gp_dual, MyMat<double> G_real, MyMat<double> &F_gp)
{	

	// F_gp -->  F^k_\gamma c_k \otimes E^\gamma
	// F_gp_dual --> F_k^\gamma c^k \otimes E_\gamma

	for (int k=0;k<dim;k++){
		for (int gamma=0;gamma<2;gamma++){

			F_gp[k][gamma] = 0.0;

			for (int p=0;p<dim;p++){
				for (int beta=0;beta<2;beta++){

					F_gp[k][gamma] = F_gp[k][gamma] + II[p][k]*G_real[beta][gamma]*F_gp_dual[p][beta];
				}
			}

		}
	}

	return 0;

}

inline PetscErrorCode Dual_Solve::get_C(MyMat<double> F_gp, MyMat <double> F_gp_dual, MyMat <double> &C)
{	
	// C_\gamma^\delta E^\gamma \otimes E_\delta

	for (int gamma=0;gamma<2;gamma++){
		for (int delta=0;delta<2;delta++){

			C[gamma][delta] = 0.0;

			for (int k=0;k<dim;k++){

				C[gamma][delta] = C[gamma][delta] + F_gp[k][gamma]*F_gp_dual[k][delta];
			}
		}
	}

	return 0;

}

inline PetscErrorCode Dual_Solve::get_L2_norm_F_diff_F_0(MyMat<double> F_guess, MyMat<double> F_guess_dual, 
											MyMat<double> F_0_real, MyMat<double> F_0_dual, double *sum_norm)
{
	*sum_norm = 0.0;

	for (int k1=0; k1<this->dim; k1++){
		for (int k2=0; k2<2; k2++){

			*sum_norm = *sum_norm + (F_guess[k1][k2] - F_0_real[k1][k2])*(F_guess_dual[k1][k2] - F_0_dual[k1][k2]);
		}
	}

	return 0;
}

inline PetscErrorCode Dual_Solve::eigenvectors_on_target_shape(MyVec <double> eig_v1, MyVec <double> eig_v2, 
														MyMat <double> F, MyMat <double> G_dual,
														MyVec <double> &Eig1 , MyVec <double> &Eig2)
{

	for (int k=0; k<2; k++){

		Eig1[k] = 0.0;
		Eig2[k] = 0.0;
	}

	for (int k=0; k<2; k++){
		for (int gamma=0; gamma<2; gamma++){
			for (int alpha=0; alpha<2; alpha++){

				Eig1[k] = Eig1[k] +  F[k][gamma]*eig_v1[alpha]*G_dual[gamma][alpha];
				Eig2[k] = Eig2[k] +  F[k][gamma]*eig_v2[alpha]*G_dual[gamma][alpha];
			}
		}
	}


	return 0;
}

inline PetscErrorCode Dual_Solve::eigenvectors_on_reference_shape(MyVec <double> eig_v1, MyVec <double> eig_v2,
															MyMat <double> E_dual_el,int i,
															MyVec <double> &Eig1, MyVec <double> &Eig2)
{
	
	for (int k=0; k<2; k++){

		Eig1[k] = 0.0;
		Eig2[k] = 0.0;

		for (int alpha=0; alpha<2; alpha++){

			Eig1[k] = Eig1[k] +  eig_v1[alpha]*E_dual_el[i][grid->X_DOF*alpha+k];
			Eig2[k] = Eig2[k] +  eig_v2[alpha]*E_dual_el[i][grid->X_DOF*alpha+k];
		}
	}

	return 0;
}



PetscErrorCode Dual_Solve::delete_petsc_objects_newton_raphson()
{

	// delete objects
	// MatDestroy(&K_xi);	VecDestroy(&F_xi); VecDestroy(&xi_vec);

	// KSPDestroy(&ksp_xi); VecScatterDestroy(&ctx_xi); VecDestroy(&xi_vec_SEQ);

	// deleting objects
	VecDestroy(&rhs_vec);	VecDestroy(&dual_vec);

	VecDestroy(&delta_dual_vec);

	VecScatterDestroy(&ctx_dual); VecDestroy(&dual_vec_SEQ);

	MatDestroy(&K_xdef); MatDestroy(&K_dual);

	VecDestroy(&F_xdef);  VecDestroy(&xdef_vec);

	MatDestroy(&Ke_local); VecDestroy(&Fe_local); VecDestroy(&local_vec); KSPDestroy(&ksp_local);

	KSPDestroy(&ksp_xdef); KSPDestroy(&ksp_dual);

	VecScatterDestroy(&ctx_xdef); VecDestroy(&xdef_vec_SEQ);

	return 0;
}


PetscErrorCode Dual_Solve::delete_petsc_objects_gradient_flow()
{

	// delete objects
	// MatDestroy(&K_xi);	VecDestroy(&F_xi); VecDestroy(&xi_vec);

	// KSPDestroy(&ksp_xi); VecScatterDestroy(&ctx_xi); VecDestroy(&xi_vec_SEQ);

	// deleting objects
	VecDestroy(&rhs_vec);	VecDestroy(&dual_vec);

	VecDestroy(&mass_vec); VecDestroy(&rhs_vec_final);

	VecScatterDestroy(&ctx_dual); VecDestroy(&dual_vec_SEQ);

	MatDestroy(&K_xdef); MatDestroy(&K_dual);

	VecDestroy(&F_xdef);  VecDestroy(&xdef_vec);

	MatDestroy(&Ke_local); VecDestroy(&Fe_local); VecDestroy(&local_vec); KSPDestroy(&ksp_local);

	KSPDestroy(&ksp_xdef); 

	VecScatterDestroy(&ctx_xdef); VecDestroy(&xdef_vec_SEQ);

	return 0;
}

