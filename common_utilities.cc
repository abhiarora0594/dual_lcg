#include "common_utilities.h"

Common_Utilities::Common_Utilities(Fem *fem, Grid *grid, Dual_Solve *dual_solve)
{	

	this->fem = fem;
	this->grid = grid;
    this->dual_solve = dual_solve; 

}

Common_Utilities::~Common_Utilities()
{

}

PetscErrorCode Common_Utilities::mysolve_local_linear_system(int n, MyMat <double> k_local,
														MyVec <double> f_local, MyVec <double> &x_local)
{

	// initialize L and U
	MyMat <double> L(n,MyVec <double>(n));
	MyMat <double> U(n,MyVec <double>(n));

	MyVec <double> y(n);

	for (int i=0;i<n;i++){
		x_local[i] = 0.0;
	}

	double sum;

	for(int i=0;i<n;i++){
	 	U[i][i] = 1.0; // for Crout reduction
	}

	// Do the LU decomposition using the Crout reduction

	for(int i=0;i<n;i++) 
	{
		// loop over pairs of L columns and U rows. The three levels of loop indicate n^3 behavior

		// first, column i of L
	 
		for(int j=i;j<n;j++)
		{ // i is the row
			
			sum=0.0;

			for(int k=0;k<=i-1;k++)
			{
				sum = sum + L[j][k]*U[k][i];
			}

			L[j][i] = k_local[j][i] - sum;
		}
		
		// second, row i of U
		   
		for(int j=i+1;j<n;j++)
		{ // j is the column
		   	
		   	sum = 0.0;

		   	for(int k=0;k<=i-1;k++)
		   	{
		    	sum = sum + L[i][k]*U[k][j];
			}

			if (std::abs(L[i][i]) < 1e-3){
				throw "LU decomposition division by zero";
			}
		    U[i][j] = (k_local[i][j] - sum)/L[i][i];
			
		}
	}

	// solve the system of equations

	// first, find the y vector

	y[0] = f_local[0]/L[0][0];
	
	for(int i=1;i<n;i++) 
	{
	    sum = 0.0;
	    
	    for(int j=0;j<i;j++)
	    {
	   		sum = sum + L[i][j]*y[j];
	    }	

	   	if (std::abs(L[i][i]) < 1e-4){
			throw "LU decomposition division by zero";
		}
	   	y[i] = (f_local[i]-sum)/L[i][i];
	
	}

	// second, find the x vector

	x_local[n-1] = y[n-1];

	int j;

	for(int i=1;i<n;i++) 
	{
		j = n-i-1;
	   
		sum = 0.0;
	   
	    for(int k=j+1;k<n;k++){
	   		sum = sum + U[j][k]*x_local[k];
	    }
	   	
	   	x_local[j] = y[j] - sum;
	}

	// for (int i=0; i<n; i++)
	// {
	// 	std::cout << "x_local is " << x_local[i] << std::endl;
	// }	

   return 0;
}


PetscErrorCode Common_Utilities::petsc_local_linear_system(int n, MyMat <double> k_local,
														MyVec <double> f_local, MyVec <double> &x_local)
{

	PetscErrorCode ierr;

	PetscReal rnorm;
	PetscInt its;

	int indc[n];
	double ke_vec[n*n];

	double l_vec[n];
	double f_vec[n];

	for (int i=0;i<n;i++){
		indc[i] = i;
		
		for (int j=0;j<n;j++){
			ke_vec[i*n+j] = k_local[i][j];
		}
		f_vec[i] = f_local[i];
	}

	ierr = VecSetValues(dual_solve->Fe_local,n,indc,f_vec,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValues(dual_solve->Ke_local,n,indc,n,indc,ke_vec,INSERT_VALUES); CHKERRQ(ierr);

	// final assembly of global vector (due to multiple proccess)
	ierr = VecAssemblyBegin(dual_solve->Fe_local); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(dual_solve->Fe_local); CHKERRQ(ierr);

	// final assembly of global matrix (due to multiple proccess)
	ierr = MatAssemblyBegin(dual_solve->Ke_local,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(dual_solve->Ke_local,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// creating a solver context
	if (dual_solve->count == 0){
		ierr = this->ksp_mumps_solver_petsc(dual_solve->ksp_local,dual_solve->Ke_local,
									dual_solve->pc_local,dual_solve->FF_local); CHKERRQ(ierr);
	}

	//Solve the linear system
	ierr = KSPSolve(dual_solve->ksp_local,dual_solve->Fe_local,dual_solve->local_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(dual_solve->ksp_local,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(dual_solve->ksp_local,&its); CHKERRQ(ierr);

	ierr = VecGetValues(dual_solve->local_vec,n,indc,l_vec); CHKERRQ(ierr);	

	for (int i=0; i<n; i++)
	{
		x_local[i] = l_vec[i];
		std::cout << "x_local is " << x_local[i] << std::endl;
	}	

   return 0;
}

PetscErrorCode Common_Utilities::vtk_write()
{

	std::ofstream myfile_vtk;
	int el_type=9; // 9 for quad
	int node_per_el= fem->node_per_el; // 4 for quad
	std::string snapshot_time= "data at T = ";

	std::string filename = "./outputs/fields_data_";

	snapshot_time = snapshot_time + std::to_string(dual_solve->t_step);

	filename = filename + std::to_string(dual_solve->step_no) + ".vtk";

	myfile_vtk.open(filename, std::ios::out);

	myfile_vtk << "# vtk DataFile Version 2.0" << "\r\n";
	myfile_vtk << snapshot_time << "\r\n";
	myfile_vtk << "ASCII" << "\r\n";
	myfile_vtk << "DATASET UNSTRUCTURED_GRID" << "\r\n";
	myfile_vtk << "POINTS " << grid->ndof << " double" << "\r\n";

	double tmp1, tmp2, tmp3, tmp4, mod;

	if (dual_solve->step_no == 0){
		for (int i=0;i<grid->ndof;++i){
			myfile_vtk << std::fixed << std::setprecision(10) << grid->X_ref[i][0] << 
						" " << grid->X_ref[i][1] << " " << grid->X_ref[i][2] << "\r\n";
		}
	}
	else {
		for (int i=0;i<grid->ndof;++i){
			myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->x_def[i][0] << 
						" " << dual_solve->x_def[i][1] << " " << 0.0 << "\r\n";
		}
	}

	myfile_vtk << "CELLS " << grid->nel << " " << grid->nel*(node_per_el+1) << "\r\n";

	for (int i=0; i<grid->nel; ++i){
		myfile_vtk << node_per_el << " " << grid->el_conn[i][0] << " " << grid->el_conn[i][1] << " " << 
							grid->el_conn[i][2] << " " << grid->el_conn[i][3] << "\r\n";
	}

	myfile_vtk << "CELL_TYPES" << " " << grid->nel << "\r\n";
	for (int i=0;i<grid->nel;++i){
		myfile_vtk << el_type << "\r\n";	
	}

	myfile_vtk << "CELL_DATA " << grid->nel << "\r\n";
	// myfile_vtk << "SCALARS " << "norm_l2_x " << "double " << 1 << "\r\n";
	// myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	// for (int ie=0;ie<grid->nel;++ie){
	// 	myfile_vtk << std::fixed << std::setprecision(10) << norm_global_x[ie] << "\r\n";
	// }

	// myfile_vtk << "SCALARS " << "norm_l2_F " << "double " << 1 << "\r\n";
	// myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	// for (int ie=0;ie<grid->nel;++ie){
	// 	myfile_vtk << std::fixed << std::setprecision(10) << norm_global_F[ie] << "\r\n";
	// }


	myfile_vtk << "SCALARS " << "dual_eig1 " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->dual_S_mu1[ie] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "dual_eig2 " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->dual_S_mu2[ie] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "dual_F_dx " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->dual_S_F_dx[ie] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "dual_H " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->dual_S_quad_H[ie] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "dual_func_total " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->dual_S_norm[ie] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "qfactor " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->qfactor[ie] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "compatibility_check " << "double " << 1 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->compatibility_check[ie] << "\r\n";
	}

	myfile_vtk << "FIELD " << "fieldData " << 4 + 6 << "\r\n";
	myfile_vtk << "x_gp " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->x_gp[ie][0][0] << 
															" " << dual_solve->x_gp[ie][0][1] << "\r\n";
	}

	myfile_vtk << "F_gp " << dual_solve->F_DOF << " " << grid->nel << " double" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->F_gp[ie][0][0] << " " 
															<< dual_solve->F_gp[ie][0][1] << " " 
															<< dual_solve->F_gp[ie][0][2] << " " 
															<< dual_solve->F_gp[ie][0][3] << "\r\n";
	}

	// myfile_vtk << "Eig1_avg " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	// for (int ie=0;ie<grid->nel;++ie){
		
	// 	tmp1 = 0.25*(dual_solve->Eig_1_gp[ie][0][0] + dual_solve->Eig_1_gp[ie][1][0] + 
	// 						dual_solve->Eig_1_gp[ie][2][0] + dual_solve->Eig_1_gp[ie][3][0]);

	// 	tmp2 = 0.25*(dual_solve->Eig_1_gp[ie][0][1] + dual_solve->Eig_1_gp[ie][1][1] + 
	// 						dual_solve->Eig_1_gp[ie][2][1] +  dual_solve->Eig_1_gp[ie][3][1]); 

	// 	mod =  std::sqrt(tmp1*tmp1 + tmp2*tmp2);

	// 	tmp1 = tmp1/mod;
	// 	tmp2 = tmp2/mod;

	// 	myfile_vtk << std::fixed << std::setprecision(10) << tmp1 << " " << tmp2 << "\r\n";
	// }

	// myfile_vtk << "Eig2_avg " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	
	// for (int ie=0;ie<grid->nel;++ie){

	// 	tmp3 = 0.25*(dual_solve->Eig_2_gp[ie][0][0] + dual_solve->Eig_2_gp[ie][1][0] + 
	// 					dual_solve->Eig_2_gp[ie][2][0] + dual_solve->Eig_2_gp[ie][3][0]);

	// 	tmp4 = 0.25*(dual_solve->Eig_2_gp[ie][0][1] + dual_solve->Eig_2_gp[ie][1][1] + 
	// 					dual_solve->Eig_2_gp[ie][2][1] + dual_solve->Eig_2_gp[ie][3][1]);

	// 	mod =  std::sqrt(tmp3*tmp3 + tmp4*tmp4);

	// 	tmp3 = tmp3/mod;
	// 	tmp4 = tmp4/mod;		 

	// 	myfile_vtk << std::fixed << std::setprecision(10) << tmp3 << " " << tmp4 << "\r\n";
	// }

	myfile_vtk << "Eig1_1 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		
		tmp1 = dual_solve->Eig_1_gp[ie][0][0];

		tmp2 = dual_solve->Eig_1_gp[ie][0][1];

		mod =  std::sqrt(tmp1*tmp1 + tmp2*tmp2);

		tmp1 = tmp1/mod;
		tmp2 = tmp2/mod;

		myfile_vtk << std::fixed << std::setprecision(10) << tmp1 << " " << tmp2 << "\r\n";
	}

	myfile_vtk << "Eig2_1 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	
	for (int ie=0;ie<grid->nel;++ie){

		tmp3 = dual_solve->Eig_2_gp[ie][0][0];

		tmp4 = dual_solve->Eig_2_gp[ie][0][1];

		mod =  std::sqrt(tmp3*tmp3 + tmp4*tmp4);

		tmp3 = tmp3/mod;
		tmp4 = tmp4/mod;		 

		myfile_vtk << std::fixed << std::setprecision(10) << tmp3 << " " << tmp4 << "\r\n";
	}

	myfile_vtk << "Eig1_2 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		
		tmp1 = dual_solve->Eig_1_gp[ie][1][0];

		tmp2 = dual_solve->Eig_1_gp[ie][1][1];

		mod =  std::sqrt(tmp1*tmp1 + tmp2*tmp2);

		tmp1 = tmp1/mod;
		tmp2 = tmp2/mod;

		myfile_vtk << std::fixed << std::setprecision(10) << tmp1 << " " << tmp2 << "\r\n";
	}

	myfile_vtk << "Eig2_2 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	
	for (int ie=0;ie<grid->nel;++ie){

		tmp3 = dual_solve->Eig_2_gp[ie][1][0];

		tmp4 = dual_solve->Eig_2_gp[ie][1][1];

		mod =  std::sqrt(tmp3*tmp3 + tmp4*tmp4);

		tmp3 = tmp3/mod;
		tmp4 = tmp4/mod;		 

		myfile_vtk << std::fixed << std::setprecision(10) << tmp3 << " " << tmp4 << "\r\n";
	}

	myfile_vtk << "Eig1_3 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		
		tmp1 = dual_solve->Eig_1_gp[ie][2][0];

		tmp2 = dual_solve->Eig_1_gp[ie][2][1];

		mod =  std::sqrt(tmp1*tmp1 + tmp2*tmp2);

		tmp1 = tmp1/mod;
		tmp2 = tmp2/mod;

		myfile_vtk << std::fixed << std::setprecision(10) << tmp1 << " " << tmp2 << "\r\n";
	}

	myfile_vtk << "Eig2_3 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	
	for (int ie=0;ie<grid->nel;++ie){

		tmp3 = dual_solve->Eig_2_gp[ie][2][0];

		tmp4 = dual_solve->Eig_2_gp[ie][2][1];

		mod =  std::sqrt(tmp3*tmp3 + tmp4*tmp4);

		tmp3 = tmp3/mod;
		tmp4 = tmp4/mod;		 

		myfile_vtk << std::fixed << std::setprecision(10) << tmp3 << " " << tmp4 << "\r\n";
	}

	myfile_vtk << "Eig1_4 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	for (int ie=0;ie<grid->nel;++ie){
		
		tmp1 = dual_solve->Eig_1_gp[ie][3][0];

		tmp2 = dual_solve->Eig_1_gp[ie][3][1];

		mod =  std::sqrt(tmp1*tmp1 + tmp2*tmp2);

		tmp1 = tmp1/mod;
		tmp2 = tmp2/mod;

		myfile_vtk << std::fixed << std::setprecision(10) << tmp1 << " " << tmp2 << "\r\n";
	}

	myfile_vtk << "Eig2_4 " << dual_solve->dim << " " << grid->nel << " double" << "\r\n";

	
	for (int ie=0;ie<grid->nel;++ie){

		tmp3 = dual_solve->Eig_2_gp[ie][3][0];

		tmp4 = dual_solve->Eig_2_gp[ie][3][1];

		mod =  std::sqrt(tmp3*tmp3 + tmp4*tmp4);

		tmp3 = tmp3/mod;
		tmp4 = tmp4/mod;		 

		myfile_vtk << std::fixed << std::setprecision(10) << tmp3 << " " << tmp4 << "\r\n";
	}

	myfile_vtk << "POINT_DATA " << grid->ndof << "\r\n";
	myfile_vtk << "SCALARS " << "beta " << "double " << 2 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int i=0;i<grid->ndof;++i){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->beta[i][0] << " " 
															<< dual_solve->beta[i][1] << " " << "\r\n";
	}

	myfile_vtk << "SCALARS " << "A_dual " << "double " << grid->A_DOF << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int i=0;i<grid->ndof;++i){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->A_dual[i][0] << " " 
															<< dual_solve->A_dual[i][1] << " " 
															<< dual_solve->A_dual[i][2] << " " 
															<< dual_solve->A_dual[i][3] << "\r\n";
	}

	myfile_vtk << "SCALARS " << "beta_rhs " << "double " << 2 << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int i=0;i<grid->ndof;++i){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->beta_rhs[i][0] << " " 
															<< dual_solve->beta_rhs[i][1] << " " << "\r\n";
	}

	myfile_vtk << "SCALARS " << "A_dual_rhs " << "double " << grid->A_DOF << "\r\n";
	myfile_vtk << "LOOKUP_TABLE " << "default" << "\r\n";

	for (int i=0;i<grid->ndof;++i){
		myfile_vtk << std::fixed << std::setprecision(10) << dual_solve->A_dual_rhs[i][0] << " " 
															<< dual_solve->A_dual_rhs[i][1] << " " 
															<< dual_solve->A_dual_rhs[i][2] << " " 
															<< dual_solve->A_dual_rhs[i][3] << "\r\n";
	}

	myfile_vtk.close();

	return 0;	
}


PetscErrorCode Common_Utilities::restart_write()
{

	std::ofstream myfile1, myfile2, myfile3, myfile4, myfile5;

	std::string filename1 = "./restart/restart_x_gp_" + std::to_string(dual_solve->step_no) + ".txt";

	myfile1.open(filename1, std::ios::out);

	for (int ie=0;ie<grid->nel;++ie){
		
		for (int gp=0;gp<fem->ngp;gp++){
			
			myfile1 << std::fixed << std::setprecision(10) << dual_solve->x_gp[ie][gp][0] << " " 
															<< dual_solve->x_gp[ie][gp][1] << "\r\n";
		}

	}

	myfile1.close();


	std::string filename2 = "./restart/restart_F_gp_" + std::to_string(dual_solve->step_no) + ".txt";

	myfile2.open(filename2, std::ios::out);

	for (int ie=0;ie<grid->nel;++ie){

		for (int gp=0;gp<fem->ngp;gp++){

			myfile2 << std::fixed << std::setprecision(10) << dual_solve->F_gp[ie][gp][0] << " " 
															<< dual_solve->F_gp[ie][gp][1] << " " 
															<< dual_solve->F_gp[ie][gp][2] << " " 
															<< dual_solve->F_gp[ie][gp][3] << "\r\n";	
		}
	}

	myfile2.close();


	std::string filename3 = "./restart/restart_t_" + std::to_string(dual_solve->step_no) + ".txt";

	myfile3.open(filename3, std::ios::out);

	myfile3 << dual_solve->step_no << " " << std::fixed << std::setprecision(10) << dual_solve->t_step << "\r\n";	

	myfile3.close();


	std::string filename4 = "./restart/restart_x_def_" + std::to_string(dual_solve->step_no) + ".txt";

	myfile4.open(filename4, std::ios::out);

	for (int i=0;i<grid->ndof;++i){
			
		myfile4 << std::fixed << std::setprecision(10) << dual_solve->x_def[i][0] << " " 
														<< dual_solve->x_def[i][1] << "\r\n";

	}

	myfile4.close();

	std::string filename5 = "./restart/restart_dual_" + std::to_string(dual_solve->step_no) + ".txt";

	myfile5.open(filename5, std::ios::out);

	for (int i=0;i<grid->ndof;++i){
			
		myfile5 << std::fixed << std::setprecision(10) << dual_solve->beta[i][0] << " " 
														<< dual_solve->beta[i][1] << " "
														<< dual_solve->A_dual[i][0] << " "
														<< dual_solve->A_dual[i][1] << " "
														<< dual_solve->A_dual[i][2] << " "
														<< dual_solve->A_dual[i][3] << "\r\n";

	}

	myfile5.close();
	

	return 0;	
}

template<class T>
inline PetscErrorCode Common_Utilities::get_max_val(MyVec<T> v, T *val)
{

	double largest_element = 0.0; 

	typename MyVec<T>::iterator it;

	for (it = v.begin(); it != v.end(); it++)
	{
	    if(*it > largest_element)
	    {
	      largest_element = *it;
	    }
	}

	*val = largest_element;

	return 0;
}


PetscErrorCode Common_Utilities::mat_create_petsc(Mat &K,int nofdof,int a,int b)
{	

	PetscErrorCode ierr;

	ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRQ(ierr); //Creates matrix
	ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,nofdof,nofdof);CHKERRQ(ierr);	//Reserves memory
	ierr = MatSetFromOptions(K);CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(K,a,NULL);CHKERRQ(ierr);	// pre-allocation for no of non-zero entries per row
	ierr = MatMPIAIJSetPreallocation(K,a,NULL,b,NULL);CHKERRQ(ierr);
	ierr = MatSetUp(K);CHKERRQ(ierr); //Finilizes matrix creation

	return 0;
}

PetscErrorCode Common_Utilities::vec_create_petsc(Vec &F,int nofdof)
{
	PetscErrorCode ierr;

	ierr = VecCreate(PETSC_COMM_WORLD,&F); CHKERRQ(ierr);
	ierr = VecSetSizes(F,PETSC_DECIDE,nofdof); CHKERRQ(ierr);
	ierr = VecSetFromOptions(F);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode Common_Utilities::ksp_mumps_solver_petsc(KSP &ksp, Mat &KG,PC &pc, Mat &F)
{	

	PetscErrorCode ierr;

	// linear solver and its various options
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,KG,KG);CHKERRQ(ierr);

	//ierr = KSPSetType(ksp,KSPPREONLY); CHKERRQ(ierr);
	// PetscInt  ival,icntl;
	// PetscReal val;
	ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);

	ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
	// ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);

	ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = PCFactorSetUpMatSolverPackage(pc); CHKERRQ(ierr); /* call MatGetFactor() to create F */
    //ierr = PCFactorGetMatrix(pc,&F); CHKERRQ(ierr);

	ierr = KSPSetTolerances(ksp,1e-16,1e-25,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr); //Finalizes set up of solver with default options
    ierr =  KSPSetUp(ksp); CHKERRQ(ierr);

	return 0;
}