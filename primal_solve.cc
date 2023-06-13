#include "primal_solve.h"
#include "common_utilities.h"

Primal_Solve::Primal_Solve(Fem *fem, Grid *grid, Dual_Solve *dual_solve)
{	

	this->fem = fem;
	this->grid = grid;
    this->dual_solve = dual_solve; 

    count = 0;
    step_no = 0;
    
    this->lambda_1 = dual_solve->lambda_1;
    this->lambda_2 = dual_solve->lambda_2;

    dof_per_node = dual_solve->F_DOF + grid->BETA_DOF + dual_solve->dim;

    dof_per_el = this->dof_per_node*fem->node_per_el;

    tdof = grid->ndof*this->dof_per_node;

    tdof_fdef = grid->ndof*dual_solve->F_DOF;

    this->Istart = (::rank)*(grid->nel/(::size)) + ((grid->nel%(::size)) < (::rank) ? (grid->nel%(::size)) : (::rank));

  	this->Iend = Istart + grid->nel/(::size) + ((grid->nel%(::size)) > (::rank));

    beta_c = MyMat <double> (grid->ndof,MyVec<double>(2)); 
	// beta_1 and beta_2 (cartesian)

    x_def = MyMat <double> (grid->ndof,MyVec<double>(2)); 
	// x_1 and x_2 (cartesian)

    F_def = MyMat <double> (grid->ndof,MyVec<double>(dual_solve->F_DOF)); 
	// F tensor in (c_k and E^gamma basis) at nodes	= {F11,F12,F21,F22}

    primal_sol = new double [this->tdof];

  	ind_set_primal = new int [this->tdof];

    f_def_sol = new double [this->tdof_fdef];

    ind_set_fdef = new int [this->tdof_fdef];

    // index set for primal 
  	for (int i=0;i<grid->ndof;i++){
  		for (int d=0;d<this->dof_per_node;d++){

  			ind_set_primal[this->dof_per_node*i+d] = this->dof_per_node*i+d;
  		}
	}

}

Primal_Solve::~Primal_Solve()
{

    delete [] primal_sol;
    delete [] ind_set_primal;

    delete [] f_def_sol;
    delete [] ind_set_fdef;

    VecDestroy(&rhs_vec);
    VecDestroy(&primal_vec);
    VecDestroy(&delta_primal_vec);

    VecDestroy(&fdef_vec);
    VecDestroy(&F_fdef);

    VecScatterDestroy(&ctx_primal); 
	VecDestroy(&primal_vec_SEQ);

    VecScatterDestroy(&ctx_fdef); 
	VecDestroy(&fdef_vec_SEQ);

    MatDestroy(&K_primal);
    KSPDestroy(&ksp_primal); 

    MatDestroy(&K_fdef);
    KSPDestroy(&ksp_fdef); 

}

PetscErrorCode Primal_Solve::set_coupling(Common_Utilities *com_uti)
{

    this->common_utilities = com_uti;

    return 0;
}

PetscErrorCode Primal_Solve::global_newton_raphson_primal()
{

	PetscErrorCode ierr;

	residual_conv = 1.0;

    double tol = 1e-10;

    // newton solve initializations
    ierr = common_utilities->vec_create_petsc(rhs_vec,this->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(primal_vec,this->tdof); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(delta_primal_vec,this->tdof); CHKERRQ(ierr);
	ierr = common_utilities->mat_create_petsc(K_primal,this->tdof,100,100); CHKERRQ(ierr);
    
    // newton solve initializations
    ierr = VecSet(rhs_vec,0.0); CHKERRQ(ierr);
	ierr = VecSet(delta_primal_vec,0.0); CHKERRQ(ierr);
    
    // deformation gradient l2 projection initializations
    ierr = common_utilities->vec_create_petsc(fdef_vec,this->tdof_fdef); CHKERRQ(ierr);
	ierr = common_utilities->vec_create_petsc(F_fdef,this->tdof_fdef); CHKERRQ(ierr);
	ierr = common_utilities->mat_create_petsc(K_fdef,this->tdof_fdef,60,60); CHKERRQ(ierr);

    // deformation gradient l2 projection initializations
    ierr = VecSet(fdef_vec,0.0); CHKERRQ(ierr);
	ierr = VecSet(F_fdef,0.0); CHKERRQ(ierr);

    // x_def initial condition
    for (int i=0; i<grid->ndof; i++){

        this->x_def[i][0] = dual_solve->x_def[i][0];
        this->x_def[i][1] = dual_solve->x_def[i][1];
    }

    // l2 projection of deformation gradient tensor on nodes
    ierr = this->deformation_gradient_l2_projection(); CHKERRQ(ierr);

    // initial condition set for the vector of primal solve
    ierr = this->initial_condition_for_primal_solve(); CHKERRQ(ierr);

    while (residual_conv > tol)
    {

        ierr = this->global_newton_solve(); CHKERRQ(ierr);

        ierr = this->global_newton_update(); CHKERRQ(ierr);

        ierr = MatSetOption(K_primal,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);

        ierr = MatZeroEntries(K_primal); CHKERRQ(ierr);

        step_no = step_no + 1;

        count = count + 1;

        residual_conv = norm_rhs;

        if (::rank == 0){

            std::cout << "**********************************" << std::endl;
            std::cout << "step_no is " << step_no << std::endl;
            std::cout << "norm_rhs is " << norm_rhs << std::fixed << std::setprecision(15) << std::endl;
            std::cout << "**********************************" << std::endl;
        }
        

    } // while loop ends here

	return 0;

}

PetscErrorCode Primal_Solve::initial_condition_for_primal_solve()
{

    PetscErrorCode ierr;

    int indc[this->dof_per_node];
	double val[this->dof_per_node];

	for (int i=0;i<grid->ndof;i++){

		indc[0] = i*this->dof_per_node;
		indc[1] = i*this->dof_per_node + 1;

		val[0] = x_def[i][0];
		val[1] = x_def[i][1];

		for (int k1=0;k1<dual_solve->dim;k1++){
			for (int k2=0;k2<2;k2++){

				val[2+2*k1+k2] = F_def[i][2*k1+k2];

				indc[2+2*k1+k2] = i*this->dof_per_node + 2 + 2*k1 + k2;
			}
		}

        beta_c[i][0] = 0.0;
		beta_c[i][1] = 0.0;

        val[6] = beta_c[i][0];
        val[7] = beta_c[i][0];

        indc[6] = i*this->dof_per_node + 6;
        indc[7] = i*this->dof_per_node + 7;

		ierr = VecSetValues(primal_vec,this->dof_per_node,indc,val,INSERT_VALUES); CHKERRQ(ierr);

		// final assembly of global vector (due to multiple proccess)
		ierr = VecAssemblyBegin(primal_vec); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(primal_vec); CHKERRQ(ierr);

	}

    return 0;
}

PetscErrorCode Primal_Solve::global_newton_update()
{
	PetscErrorCode ierr;

	double scal = 1.0;

	ierr = VecNorm(rhs_vec, NORM_INFINITY, &norm_rhs); CHKERRQ(ierr);

	ierr = VecAXPY(primal_vec,scal,delta_primal_vec); CHKERRQ(ierr);

	ierr = this->copying_petsc_primal_vector_to_stdvector(); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode Primal_Solve::global_newton_solve()
{

    PetscErrorCode ierr;

	PetscReal rnorm;
	PetscInt its;

	MyMat<double> beta_el(fem->node_per_el,MyVec<double> (grid->BETA_DOF));

	MyMat<double> F_el(fem->node_per_el, MyVec<double>(dual_solve->F_DOF));

    MyMat<double> F_gp(fem->ngp, MyVec<double>(dual_solve->F_DOF));

    MyMat<double> xdef_el(fem->node_per_el, MyVec<double>(dual_solve->dim));

	int node_el;

	double rhs_el[this->dof_per_el]; // element rhs vector

	double **ke; // element stiffness matrix

	ke = new double* [this->dof_per_el];
	for (int i=0;i<this->dof_per_el;i++)
		ke[i]= new double [this->dof_per_el];

	double ke_vec[this->dof_per_el*this->dof_per_el];

	int indc[this->dof_per_el];

	MyVec <double> detJ_e(fem->ngp);

	MyMat <double> E_dual_el(fem->ngp,MyVec<double>(2*grid->X_DOF));

	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initializing rhs_el to be zero (passed as pointers: should not sum)
		for (int i=0;i<this->dof_per_el;i++){
			for (int j=0;j<this->dof_per_el;j++){
				ke[i][j] = 0.0;
			}
			rhs_el[i] = 0.0;
		}

		// getting fields for each elem
		for (int i=0;i<fem->node_per_el;i++){

			node_el = grid->el_conn[ie][i];
			
			for (int k=0;k<grid->BETA_DOF;k++){
				beta_el[i][k] = beta_c[node_el][k];
			}

			for (int k=0;k<dual_solve->F_DOF;k++){
				F_el[i][k] = F_def[node_el][k];
			}

			for (int k=0;k<dual_solve->dim;k++){
				xdef_el[i][k] = x_def[node_el][k];
			}
		}

        if (this->count == 0){
            
            // getting F at gps from xdef (initial_guess)
            for(int i=0;i<fem->ngp; i++){

                F_gp[i][0] = 0.0; F_gp[i][1] = 0.0; F_gp[i][2] = 0.0; F_gp[i][3] = 0.0;

                for (int p=0;p<fem->node_per_el;p++){
                    
                    F_gp[i][0] = F_gp[i][0] + fem->dpsi[i][p][0]*xdef_el[p][0];
                    F_gp[i][1] = F_gp[i][1] + fem->dpsi[i][p][1]*xdef_el[p][0];
                    F_gp[i][2] = F_gp[i][2] + fem->dpsi[i][p][0]*xdef_el[p][1];
                    F_gp[i][3] = F_gp[i][3] + fem->dpsi[i][p][1]*xdef_el[p][1];
                }
            }

        }
        else {

            // interpolating on gps from F_el
            for (int i=0; i<fem->ngp;i++){
                for (int k=0; k<dual_solve->F_DOF; k++){

                    F_gp[i][k] = 0.0;

                    for (int p=0;p<fem->node_per_el;p++){

                        F_gp[i][k] = F_gp[i][k] + fem->psi[i][p]*F_el[p][k];
                    }
                }
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

		// getting indicies for assembly of elem vector to global vector
		for (int k1=0;k1<fem->node_per_el;k1++){
            for (int k2=0; k2<this->dof_per_node;k2++){

                indc[k1*this->dof_per_node + k2] = this->dof_per_node*grid->el_conn[ie][k1] + k2;

            }
		}

		// stiff matrix and rhs vector per elem
		ierr = this->elem_stiff_force_vector_global_newton(beta_el,F_gp,E_dual_el,detJ_e,
                                                            xdef_el,rhs_el,ke); CHKERRQ(ierr);


		// adding values to global matrix and vector
		for (int i=0;i<this->dof_per_el;i++){
			for (int j=0;j<this->dof_per_el;j++){

				ke_vec[i*this->dof_per_el+j] = ke[i][j];
			}
		}	

		// assembly to global rhs vector and stiff matrix
		ierr = VecSetValues(rhs_vec,this->dof_per_el,indc,rhs_el,ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(K_primal,this->dof_per_el,indc,this->dof_per_el,indc,ke_vec,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector (due to multiple proccess)
	ierr = VecAssemblyBegin(rhs_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs_vec); CHKERRQ(ierr);

	// final assembly of global matrix (due to multiple proccess)
	ierr = MatAssemblyBegin(K_primal,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K_primal,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// boundary conditions for newton solve
	// ierr = this->global_newton_solve_bcs(K_primal,rhs_vec,delta_primal_vec); CHKERRQ(ierr);

	// creating the solver context
	if (this->count == 0){
		ierr = common_utilities->ksp_mumps_solver_petsc(ksp_primal,K_primal,pc_primal,FF_primal); CHKERRQ(ierr);
	}
		
	//Solve the linear system
	ierr = KSPSolve(ksp_primal,rhs_vec,delta_primal_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(ksp_primal,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(ksp_primal,&its); CHKERRQ(ierr);

	if (::rank == 0){
		std::cout << "primal global newton solve residual for step no " << this->step_no << " is " << rnorm << std::endl;
	}

	// delete objects and matrices
	for (int i = 0; i<this->dof_per_el; i++){
		delete[] ke[i];
    }
	delete[] ke;

    return 0;
}

PetscErrorCode Primal_Solve::copying_petsc_primal_vector_to_stdvector()
{

    PetscErrorCode ierr;

    if (this->count == 0){
		ierr = VecScatterCreateToAll(primal_vec,&ctx_primal,&primal_vec_SEQ); CHKERRQ(ierr);
	}

	ierr = VecScatterBegin(ctx_primal,primal_vec,primal_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_primal,primal_vec,primal_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(primal_vec_SEQ,this->tdof,ind_set_primal,primal_sol); CHKERRQ(ierr);

	for (int i=0;i<grid->ndof;i++){

		x_def[i][0] = primal_sol[this->dof_per_node*i];
		x_def[i][1] = primal_sol[this->dof_per_node*i+1];

		F_def[i][0] = primal_sol[this->dof_per_node*i+2];
		F_def[i][1] = primal_sol[this->dof_per_node*i+3];
		F_def[i][2] = primal_sol[this->dof_per_node*i+4];
		F_def[i][3] = primal_sol[this->dof_per_node*i+5];

        beta_c[i][0] = primal_sol[this->dof_per_node*i+6];
        beta_c[i][1] = primal_sol[this->dof_per_node*i+7];

	}// loop over nodes

    return 0;
}

PetscErrorCode Primal_Solve::elem_stiff_force_vector_global_newton(MyMat<double> beta_el, MyMat<double>F_gp,
                                                            MyMat<double>E_dual_el, MyVec<double>detJ_e,
                                                            MyMat <double> xdef_el,
                                                            double *rhs_el, double **ke)
{

    PetscErrorCode ierr;

    double fe[this->dof_per_el];

	double K[this->dof_per_el][this->dof_per_el];

    MyMat<double> F(dual_solve->dim,MyVec<double>(2));

	MyMat<double> F_dual(dual_solve->dim,MyVec<double>(2));

	MyMat<double> G_dual(2,MyVec<double>(2));

	MyMat<double> G_real(2,MyVec<double>(2));

	MyMat<double> C(2,MyVec<double>(2));

	double tr_C, det_C;

	double mu1, mu2;

    int tmp;

	MyVec <double> eig_v1(2);

	MyVec <double> eig_v2(2);

	MyVec <double> d_mu1_F(dual_solve->F_DOF);

	MyVec <double> d_mu2_F(dual_solve->F_DOF);

	MyMat <double> dd_mu1_FF(dual_solve->F_DOF,MyVec <double> (dual_solve->F_DOF));

	MyMat <double> dd_mu2_FF(dual_solve->F_DOF,MyVec <double> (dual_solve->F_DOF));

    double beta_1_gp, beta_2_gp; 

    double dx1_zeta_1, dx1_zeta_2;

    double dx2_zeta_1, dx2_zeta_2;

    double scal1, scal2, scal3, scal4, scal5, scal6, scal7, scal8;

    double scal9, scal10, scal11, scal12;

    MyMat <double> dF_zeta (dual_solve->F_DOF, MyVec <double> (2));

    MyMat <double> dE1_dual_zeta(grid->X_DOF, MyVec<double> (2));
    MyMat <double> dE2_dual_zeta(grid->X_DOF, MyVec<double> (2));

    MyVec <double> F11_res1(fem->node_per_el);
    MyVec <double> F12_res1(fem->node_per_el);
    MyVec <double> F21_res1(fem->node_per_el);
    MyVec <double> F22_res1(fem->node_per_el);

    MyVec <double> F11_res2(fem->node_per_el);
    MyVec <double> F12_res2(fem->node_per_el);
    MyVec <double> F21_res2(fem->node_per_el);
    MyVec <double> F22_res2(fem->node_per_el);

    MyVec <double> F11_res3(fem->node_per_el);
    MyVec <double> F12_res3(fem->node_per_el);
    MyVec <double> F21_res3(fem->node_per_el);
    MyVec <double> F22_res3(fem->node_per_el);

    MyVec <double> F11_res4(fem->node_per_el);
    MyVec <double> F12_res4(fem->node_per_el);
    MyVec <double> F21_res4(fem->node_per_el);
    MyVec <double> F22_res4(fem->node_per_el);

    MyVec <double> v1(grid->X_DOF);
    MyVec <double> v2(grid->X_DOF);
    MyVec <double> v3(grid->X_DOF);
    MyVec <double> v4(grid->X_DOF);

    MyVec <double> E_dual(grid->X_DOF);
    MyVec <double> E_dual2(grid->X_DOF);
    MyVec <double> E_dual3(grid->X_DOF);

    MyVec <double> F11_stiff1(dual_solve->F_DOF);
    MyVec <double> F12_stiff1(dual_solve->F_DOF);
    MyVec <double> F21_stiff1(dual_solve->F_DOF);
    MyVec <double> F22_stiff1(dual_solve->F_DOF);

    MyMat <double> F11_stiff2(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));
    MyMat <double> F12_stiff2(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));
    MyMat <double> F21_stiff2(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));
    MyMat <double> F22_stiff2(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));

    MyMat <double> F11_stiff3(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));
    MyMat <double> F12_stiff3(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));
    MyMat <double> F21_stiff3(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));
    MyMat <double> F22_stiff3(fem->node_per_el, MyVec <double> (dual_solve->F_DOF));

    // derivative of F_gp and E_dual wrt \xi coordinates at gpts
    for (int j=0;j<dual_solve->F_DOF; j++){
        for (int k=0;k<fem->ngp;k++){
            dF_zeta[j][0] = dF_zeta[j][0] + fem->dpsi_c[k][0]*F_gp[k][j];
            dF_zeta[j][1] = dF_zeta[j][1] + fem->dpsi_c[k][1]*F_gp[k][j]; 
        }
    }

    for (int j=0;j<grid->X_DOF;j++){
        for (int k=0;k<fem->ngp;k++){
            dE1_dual_zeta[j][0] = dE1_dual_zeta[j][0] + fem->dpsi_c[k][0]*E_dual_el[k][j];
            dE1_dual_zeta[j][1] = dE1_dual_zeta[j][1] + fem->dpsi_c[k][1]*E_dual_el[k][j];

            dE2_dual_zeta[j][0] = dE2_dual_zeta[j][0] + fem->dpsi_c[k][0]*E_dual_el[k][grid->X_DOF + j];
            dE2_dual_zeta[j][1] = dE2_dual_zeta[j][1] + fem->dpsi_c[k][1]*E_dual_el[k][grid->X_DOF + j];
        }
    }
    

    // loop over the gauss points begin
    for (int i=0;i<fem->ngp;i++)
    {

        for (int p=0;p<this->dof_per_el;p++){
			for (int q=0;q<this->dof_per_el;q++){
                K[p][q] = 0.0; // zeroing before summing
            }
            fe[p] = 0.0; // zeroing before summing
        }

        // zeroing values before summing
        beta_1_gp = 0.0;
        beta_2_gp = 0.0;
        dx1_zeta_1 = 0.0;
        dx1_zeta_2 = 0.0;
        dx2_zeta_1 = 0.0;
        dx2_zeta_2 = 0.0;

        for (int k1=0;k1<dual_solve->dim;k1++){
			for(int k2=0;k2<2;k2++){

				tmp = k1*2 + k2;
				F[k1][k2] = F_gp[i][tmp];
			}
		}

        ierr = dual_solve->get_G_dual(i,E_dual_el,G_dual); CHKERRQ(ierr);

		ierr = fem->my_mat_inverse(G_dual,G_real); CHKERRQ(ierr);		

		ierr = dual_solve->tranform_F_to_dual(F,G_dual,F_dual); CHKERRQ(ierr);

		ierr = dual_solve->get_C(F,F_dual,C); CHKERRQ(ierr);

		ierr = dual_solve->trace_of_C(C,&tr_C); CHKERRQ(ierr);

		ierr = dual_solve->determinant_of_C(F,G_real,&det_C); CHKERRQ(ierr);
		
		ierr = dual_solve->get_invariants(tr_C,det_C,&mu1,&mu2); CHKERRQ(ierr);

		ierr = dual_solve->get_eigenvectors_of_C(C,mu1,mu2,G_dual,eig_v1,eig_v2); CHKERRQ(ierr);

		ierr = dual_solve->first_derivative_of_invariants(eig_v1,eig_v2,G_dual,F,d_mu1_F,d_mu2_F); CHKERRQ(ierr);

		ierr = dual_solve->second_derivative_of_invariants(eig_v1,eig_v2,mu1,mu2,tr_C,C,G_dual,G_real,
                                                            F,dd_mu1_FF,dd_mu2_FF); CHKERRQ(ierr);

        // constraints evaluation at gpts    
        for (int k=0;k<fem->node_per_el;k++){
            beta_1_gp = beta_1_gp + beta_el[k][0]*fem->psi[i][k];
            beta_2_gp = beta_2_gp + beta_el[k][1]*fem->psi[i][k];
        }

        // derivative of \bf{x} wrt \xi coordinates at gpts
        for (int k=0;k<fem->node_per_el;k++){
            dx1_zeta_1 = dx1_zeta_1 + fem->dpsi[i][k][0]*xdef_el[k][0];
            dx1_zeta_2 = dx1_zeta_2 + fem->dpsi[i][k][1]*xdef_el[k][0];

            dx2_zeta_1 = dx2_zeta_1 + fem->dpsi[i][k][0]*xdef_el[k][1];
            dx2_zeta_2 = dx2_zeta_2 + fem->dpsi[i][k][1]*xdef_el[k][1];
        }


        // calculation of some residual terms
        for (int p=0;p<fem->node_per_el; p++)
        {   
            F11_res1[p] = 0.0;
            F12_res1[p] = 0.0;
            F21_res1[p] = 0.0;
            F22_res1[p] = 0.0;

            F11_res2[p] = 0.0;
            F12_res2[p] = 0.0;
            F21_res2[p] = 0.0;
            F22_res2[p] = 0.0;

            F11_res3[p] = 0.0;
            F12_res3[p] = 0.0;
            F21_res3[p] = 0.0;
            F22_res3[p] = 0.0;

            F11_res4[p] = 0.0;
            F12_res4[p] = 0.0;
            F21_res4[p] = 0.0;
            F22_res4[p] = 0.0;

            for (int delta=0; delta<2; delta++){
                for (int alpha=0; alpha<2; alpha++){
                    for (int beta=0; beta<2; beta++){
                        
                        // 1st term in residual from grad F
                        F11_res1[p] = F11_res1[p] + dF_zeta[delta][alpha]*fem->dpsi[i][p][beta]*
                                                            G_dual[delta][0]*G_dual[alpha][beta];
                        
                        F12_res1[p] = F12_res1[p] + dF_zeta[delta][alpha]*fem->dpsi[i][p][beta]*
                                                            G_dual[delta][1]*G_dual[alpha][beta];

                        F21_res1[p] = F21_res1[p] + dF_zeta[dual_solve->dim + delta][alpha]*fem->dpsi[i][p][beta]*
                                                            G_dual[delta][0]*G_dual[alpha][beta];

                        F22_res1[p] = F22_res1[p] + dF_zeta[dual_solve->dim + delta][alpha]*fem->dpsi[i][p][beta]*
                                                            G_dual[delta][1]*G_dual[alpha][beta];

                        if (delta == 0)
                        {
                            v1[0] = dE1_dual_zeta[0][alpha];
                            v1[1] = dE1_dual_zeta[1][alpha];
                            v1[2] = dE1_dual_zeta[2][alpha];

                        }
                        else if (delta == 1)
                        {
                            v1[0] = dE2_dual_zeta[0][alpha];
                            v1[1] = dE2_dual_zeta[1][alpha];
                            v1[2] = dE2_dual_zeta[2][alpha];
                        }

                        E_dual[0] = E_dual_el[i][grid->X_DOF*delta];
                        E_dual[1] = E_dual_el[i][grid->X_DOF*delta+1];
                        E_dual[2] = E_dual_el[i][grid->X_DOF*delta+2];

                        v2[0] = dE1_dual_zeta[0][beta];
                        v2[1] = dE1_dual_zeta[1][beta];
                        v2[2] = dE1_dual_zeta[2][beta];

                        v3[0] = dE2_dual_zeta[0][beta];
                        v3[1] = dE2_dual_zeta[1][beta];
                        v3[2] = dE2_dual_zeta[2][beta];

                        scal1 = 0.0;
                        scal2 = 0.0;
                        scal3 = 0.0;
                        scal4 = 0.0;
                        scal5 = 0.0;
                        scal6 = 0.0;

                        E_dual2[0] = E_dual_el[i][0];
                        E_dual2[1] = E_dual_el[i][1];
                        E_dual2[2] = E_dual_el[i][2];

                        E_dual3[0] = E_dual_el[i][grid->X_DOF];
                        E_dual3[1] = E_dual_el[i][grid->X_DOF+1];
                        E_dual3[2] = E_dual_el[i][grid->X_DOF+2];

                        ierr = fem->dot_product(v1,v2,&scal1); CHKERRQ(ierr);
                        ierr = fem->dot_product(v1,v3,&scal2); CHKERRQ(ierr);

                        ierr = fem->dot_product(E_dual,v2,&scal3); CHKERRQ(ierr);
                        ierr = fem->dot_product(E_dual,v3,&scal4); CHKERRQ(ierr);

                        ierr = fem->dot_product(v1,E_dual2,&scal5); CHKERRQ(ierr);
                        ierr = fem->dot_product(v1,E_dual3,&scal6); CHKERRQ(ierr);

                        // 2nd term in residual from grad F
                        F11_res2[p] = F11_res2[p] + fem->psi[i][p]*F[0][delta]*G_dual[alpha][beta]*scal1;
                        F12_res2[p] = F12_res2[p] + fem->psi[i][p]*F[0][delta]*G_dual[alpha][beta]*scal2;
                        F21_res2[p] = F21_res2[p] + fem->psi[i][p]*F[1][delta]*G_dual[alpha][beta]*scal1;
                        F22_res2[p] = F22_res2[p] + fem->psi[i][p]*F[1][delta]*G_dual[alpha][beta]*scal2;

                        // 3rd term in residual from grad F
                        F11_res3[p] = F11_res3[p] + dF_zeta[delta][alpha]*fem->psi[i][p]*
                                                            G_dual[alpha][beta]*scal3;
                        
                        F12_res3[p] = F12_res3[p] + dF_zeta[delta][alpha]*fem->psi[i][p]*
                                                            G_dual[alpha][beta]*scal4;
                        
                        F21_res3[p] = F21_res3[p] + dF_zeta[dual_solve->dim + delta][alpha]*fem->psi[i][p]*
                                                            G_dual[alpha][beta]*scal3;
                        
                        F22_res3[p] = F22_res3[p] + dF_zeta[dual_solve->dim + delta][alpha]*fem->psi[i][p]*
                                                            G_dual[alpha][beta]*scal4;

                        // 4th term in residual from grad F
                        F11_res4[p] = F11_res4[p] + fem->dpsi[i][p][beta]*F[0][delta]*G_dual[alpha][beta]*scal5;
                        F12_res4[p] = F12_res4[p] + fem->dpsi[i][p][beta]*F[0][delta]*G_dual[alpha][beta]*scal6;
                        F21_res4[p] = F21_res4[p] + fem->dpsi[i][p][beta]*F[1][delta]*G_dual[alpha][beta]*scal5;
                        F22_res4[p] = F22_res4[p] + fem->dpsi[i][p][beta]*F[1][delta]*G_dual[alpha][beta]*scal6;

                    }
                }
            }
        }


        // calculation of some coeffs for stiffness matrix in elliptic term
        for (int k=0; k<dual_solve->F_DOF; k++){   
            F11_stiff1[k] = 0.0;
            F12_stiff1[k] = 0.0;
            F21_stiff1[k] = 0.0;
            F22_stiff1[k] = 0.0;

            for (int p=0; p<fem->node_per_el; p++){
                F11_stiff2[p][k] = 0.0;
                F12_stiff2[p][k] = 0.0;
                F21_stiff2[p][k] = 0.0;
                F22_stiff2[p][k] = 0.0;

                F11_stiff3[p][k] = 0.0;
                F12_stiff3[p][k] = 0.0;
                F21_stiff3[p][k] = 0.0;
                F22_stiff3[p][k] = 0.0;
            }
        }    

        for (int alpha=0; alpha<2; alpha++){
            for (int beta=0; beta<2; beta++){

                scal1 = 0.0; scal2 = 0.0; scal3 = 0.0; scal4 = 0.0;

                scal5 = 0.0; scal6 = 0.0; scal7 = 0.0; scal8 = 0.0;

                scal9 = 0.0; scal10 = 0.0; scal11 = 0.0; scal12 = 0.0;

                v1[0] = dE1_dual_zeta[0][alpha];
                v1[1] = dE1_dual_zeta[1][alpha];
                v1[2] = dE1_dual_zeta[2][alpha];

                v2[0] = dE1_dual_zeta[0][beta];
                v2[1] = dE1_dual_zeta[1][beta];
                v2[2] = dE1_dual_zeta[2][beta];

                v3[0] = dE2_dual_zeta[0][alpha];
                v3[1] = dE2_dual_zeta[1][alpha];
                v3[2] = dE2_dual_zeta[2][alpha];

                v4[0] = dE2_dual_zeta[0][beta];
                v4[1] = dE2_dual_zeta[1][beta];
                v4[2] = dE2_dual_zeta[2][beta];

                E_dual2[0] = E_dual_el[i][0];
                E_dual2[1] = E_dual_el[i][1];
                E_dual2[2] = E_dual_el[i][2];

                E_dual3[0] = E_dual_el[i][grid->X_DOF];
                E_dual3[1] = E_dual_el[i][grid->X_DOF+1];
                E_dual3[2] = E_dual_el[i][grid->X_DOF+2];

                ierr = fem->dot_product(v1,v2,&scal1); CHKERRQ(ierr);
                ierr = fem->dot_product(v3,v2,&scal2); CHKERRQ(ierr);

                ierr = fem->dot_product(v1,v4,&scal3); CHKERRQ(ierr);
                ierr = fem->dot_product(v3,v4,&scal4); CHKERRQ(ierr);

                ierr = fem->dot_product(E_dual2,v2,&scal5); CHKERRQ(ierr);
                ierr = fem->dot_product(E_dual3,v2,&scal6); CHKERRQ(ierr);

                ierr = fem->dot_product(E_dual2,v4,&scal7); CHKERRQ(ierr);
                ierr = fem->dot_product(E_dual3,v4,&scal8); CHKERRQ(ierr);

                ierr = fem->dot_product(E_dual2,v1,&scal9); CHKERRQ(ierr);
                ierr = fem->dot_product(E_dual3,v1,&scal10); CHKERRQ(ierr);

                ierr = fem->dot_product(E_dual2,v3,&scal11); CHKERRQ(ierr);
                ierr = fem->dot_product(E_dual3,v3,&scal12); CHKERRQ(ierr);

                // term 1 
                // F11
                F11_stiff1[0] = F11_stiff1[0] + scal1*G_dual[alpha][beta];
                F11_stiff1[1] = F11_stiff1[1] + scal2*G_dual[alpha][beta];

                // F12
                F12_stiff1[0] = F12_stiff1[0] + scal3*G_dual[alpha][beta];
                F12_stiff1[1] = F12_stiff1[1] + scal4*G_dual[alpha][beta];

                // F21
                F21_stiff1[2] = F21_stiff1[2] + scal1*G_dual[alpha][beta];
                F21_stiff1[3] = F21_stiff1[3] + scal2*G_dual[alpha][beta];

                // F22
                F22_stiff1[2] = F22_stiff1[2] + scal3*G_dual[alpha][beta];
                F22_stiff1[3] = F22_stiff1[3] + scal4*G_dual[alpha][beta];

                // term 2
                for (int p=0; p<fem->node_per_el; p++){
                    F11_stiff2[p][0] = F11_stiff2[p][0] + scal5*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];
                    F11_stiff2[p][1] = F11_stiff2[p][1] + scal6*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];

                    F12_stiff2[p][0] = F12_stiff2[p][0] + scal7*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];
                    F12_stiff2[p][1] = F12_stiff2[p][1] + scal8*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];

                    F21_stiff2[p][2] = F21_stiff2[p][2] + scal5*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];
                    F21_stiff2[p][3] = F21_stiff2[p][3] + scal6*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];

                    F22_stiff2[p][2] = F22_stiff2[p][2] + scal7*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];
                    F22_stiff2[p][3] = F22_stiff2[p][3] + scal8*G_dual[alpha][beta]*fem->dpsi[i][p][alpha];
                }

                 // term 3
                for (int p=0; p<fem->node_per_el; p++){

                    F11_stiff3[p][0] = F11_stiff3[p][0] + scal9*G_dual[alpha][beta]*fem->dpsi[i][p][beta];
                    F11_stiff3[p][1] = F11_stiff3[p][1] + scal11*G_dual[alpha][beta]*fem->dpsi[i][p][beta];

                    F12_stiff3[p][0] = F12_stiff3[p][0] + scal10*G_dual[alpha][beta]*fem->dpsi[i][p][beta];
                    F12_stiff3[p][1] = F12_stiff3[p][1] + scal12*G_dual[alpha][beta]*fem->dpsi[i][p][beta];

                    F21_stiff3[p][2] = F21_stiff3[p][2] + scal9*G_dual[alpha][beta]*fem->dpsi[i][p][beta];
                    F21_stiff3[p][3] = F21_stiff3[p][3] + scal11*G_dual[alpha][beta]*fem->dpsi[i][p][beta];

                    F22_stiff3[p][2] = F22_stiff3[p][2] + scal10*G_dual[alpha][beta]*fem->dpsi[i][p][beta];
                    F22_stiff3[p][3] = F22_stiff3[p][3] + scal12*G_dual[alpha][beta]*fem->dpsi[i][p][beta];
                }
            }
        }


        for (int p=0; p<fem->node_per_el; p++)
        {
            for (int q=0; q<fem->node_per_el; q++)
            {
                //stiffness matrix
                //*********** x_1 ***********************
                K[p*this->dof_per_node][q*this->dof_per_node] = beta_1_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[0][0]
                                                                        + fem->dpsi[i][q][1]*dd_mu1_FF[1][0]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[0][1] + 
                                                                        fem->dpsi[i][q][1]*dd_mu1_FF[1][1])) +
                                                                    beta_2_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[0][0]
                                                                        + fem->dpsi[i][q][1]*dd_mu2_FF[1][0]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[0][1] + 
                                                                        fem->dpsi[i][q][1]*dd_mu2_FF[1][1])) +
                                                                        
                                                                    PEN*(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] +
                                                                        fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                        fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                        fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);
                                                                        
                K[p*this->dof_per_node][q*this->dof_per_node+1] = beta_1_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[2][0]
                                                                        + fem->dpsi[i][q][1]*dd_mu1_FF[3][0]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[2][1] + 
                                                                        fem->dpsi[i][q][1]*dd_mu1_FF[3][1])) +
                                                                    beta_2_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[2][0]
                                                                        + fem->dpsi[i][q][1]*dd_mu2_FF[3][0]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[2][1] + 
                                                                        fem->dpsi[i][q][1]*dd_mu2_FF[3][1]));

                K[p*this->dof_per_node][q*this->dof_per_node+2] = -PEN*fem->psi[i][q]*(fem->dpsi[i][p][0]*G_dual[0][0] + 
                                                                            fem->dpsi[i][p][1]*G_dual[0][1]);
                K[p*this->dof_per_node][q*this->dof_per_node+3] = -PEN*fem->psi[i][q]*(fem->dpsi[i][p][0]*G_dual[1][0] + 
                                                                            fem->dpsi[i][p][1]*G_dual[1][1]);
                K[p*this->dof_per_node][q*this->dof_per_node+4] = 0.0;
                K[p*this->dof_per_node][q*this->dof_per_node+5] = 0.0;

                K[p*this->dof_per_node][q*this->dof_per_node+6] = (fem->dpsi[i][p][0]*d_mu1_F[0] + 
                                                                        fem->dpsi[i][p][1]*d_mu1_F[1])*fem->psi[i][q];

                K[p*this->dof_per_node][q*this->dof_per_node+7] = (fem->dpsi[i][p][0]*d_mu2_F[0] + 
                                                                        fem->dpsi[i][p][1]*d_mu2_F[1])*fem->psi[i][q];


                // ***************** x_2 *************************
                K[p*this->dof_per_node+1][q*this->dof_per_node] = beta_1_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[0][2]
                                                                        + fem->dpsi[i][q][1]*dd_mu1_FF[1][2]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[0][3] + 
                                                                        fem->dpsi[i][q][1]*dd_mu1_FF[1][3])) +
                                                                    beta_2_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[0][2]
                                                                        + fem->dpsi[i][q][1]*dd_mu2_FF[1][2]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[0][3] + 
                                                                        fem->dpsi[i][q][1]*dd_mu2_FF[1][3]));

                K[p*this->dof_per_node+1][q*this->dof_per_node+1] = beta_1_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[2][2]
                                                                        + fem->dpsi[i][q][1]*dd_mu1_FF[3][2]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu1_FF[2][3] + 
                                                                        fem->dpsi[i][q][1]*dd_mu1_FF[3][3])) +
                                                                    beta_2_gp*(fem->dpsi[i][p][0]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[2][2]
                                                                        + fem->dpsi[i][q][1]*dd_mu2_FF[3][2]) + 
                                                                        fem->dpsi[i][p][1]*
                                                                        (fem->dpsi[i][q][0]*dd_mu2_FF[2][3] + 
                                                                        fem->dpsi[i][q][1]*dd_mu2_FF[3][3])) +
                                                                        
                                                                    PEN*(fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] +
                                                                        fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                        fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                        fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+1][q*this->dof_per_node+2] = 0.0;
                K[p*this->dof_per_node+1][q*this->dof_per_node+3] = 0.0;
                K[p*this->dof_per_node+1][q*this->dof_per_node+4] = -PEN*fem->psi[i][q]*(fem->dpsi[i][p][0]*G_dual[0][0] + 
                                                                            fem->dpsi[i][p][1]*G_dual[0][1]);

                K[p*this->dof_per_node+1][q*this->dof_per_node+5] = -PEN*fem->psi[i][q]*(fem->dpsi[i][p][0]*G_dual[1][0] + 
                                                                            fem->dpsi[i][p][1]*G_dual[1][1]);

                K[p*this->dof_per_node+1][q*this->dof_per_node+6] = (fem->dpsi[i][p][0]*d_mu1_F[2] + 
                                                                        fem->dpsi[i][p][1]*d_mu1_F[3])*fem->psi[i][q];

                K[p*this->dof_per_node+1][q*this->dof_per_node+7] = (fem->dpsi[i][p][0]*d_mu2_F[2] + 
                                                                        fem->dpsi[i][p][1]*d_mu2_F[3])*fem->psi[i][q];

                // ******************** F_11 ***********************
                K[p*this->dof_per_node+2][q*this->dof_per_node] = -PEN*fem->psi[i][p]*(fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                                fem->dpsi[i][q][1]*G_dual[1][0]);

                K[p*this->dof_per_node+2][q*this->dof_per_node+1] = 0.0;

                K[p*this->dof_per_node+2][q*this->dof_per_node+2] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[0][0]) + 
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F11_stiff1[0] +
                                                                        EPS*fem->psi[i][p]*F11_stiff2[q][0] +  
                                                                        EPS*fem->psi[i][q]*F11_stiff3[p][0] + 
                                                                        EPS*G_dual[0][0]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+2][q*this->dof_per_node+3] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[1][0]) +
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F11_stiff1[1] +
                                                                        EPS*fem->psi[i][p]*F11_stiff2[q][1] + 
                                                                        EPS*fem->psi[i][q]*F11_stiff3[p][1] + 
                                                                        EPS*G_dual[1][0]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+2][q*this->dof_per_node+4] = 0.0;
                K[p*this->dof_per_node+2][q*this->dof_per_node+5] = 0.0;
                K[p*this->dof_per_node+2][q*this->dof_per_node+6] = 0.0;
                K[p*this->dof_per_node+2][q*this->dof_per_node+7] = 0.0;

                //********************** F_12 ***********************
                K[p*this->dof_per_node+3][q*this->dof_per_node] = -PEN*fem->psi[i][p]*(fem->dpsi[i][q][0]*G_dual[0][1] + 
                                                                                fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+3][q*this->dof_per_node+1] = 0.0;

                K[p*this->dof_per_node+3][q*this->dof_per_node+2] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[0][1]) + 
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F12_stiff1[0] +
                                                                        EPS*fem->psi[i][p]*F12_stiff2[q][0] + 
                                                                        EPS*fem->psi[i][q]*F12_stiff3[q][0] +  
                                                                         EPS*G_dual[0][1]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+3][q*this->dof_per_node+3] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[1][1]) + 
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F12_stiff1[1] + 
                                                                        EPS*fem->psi[i][p]*F12_stiff2[q][1] +
                                                                        EPS*fem->psi[i][q]*F12_stiff3[q][1] + 
                                                                         EPS*G_dual[1][1]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+3][q*this->dof_per_node+4] = 0.0;
                K[p*this->dof_per_node+3][q*this->dof_per_node+5] = 0.0;
                K[p*this->dof_per_node+3][q*this->dof_per_node+6] = 0.0;
                K[p*this->dof_per_node+3][q*this->dof_per_node+7] = 0.0;

                //********************** F_21 ***********************
                K[p*this->dof_per_node+4][q*this->dof_per_node] = 0.0;
                K[p*this->dof_per_node+4][q*this->dof_per_node+1] = -PEN*fem->psi[i][p]*(fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                                fem->dpsi[i][q][1]*G_dual[1][0]);
                K[p*this->dof_per_node+4][q*this->dof_per_node+2] = 0.0;
                K[p*this->dof_per_node+4][q*this->dof_per_node+3] = 0.0;
                K[p*this->dof_per_node+4][q*this->dof_per_node+4] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[0][0]) +
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F21_stiff1[2] +
                                                                        EPS*fem->psi[i][p]*F21_stiff2[q][2] + 
                                                                        EPS*fem->psi[i][q]*F21_stiff3[p][2] + 
                                                                         EPS*G_dual[0][0]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+4][q*this->dof_per_node+5] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[1][0]) + 
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F21_stiff1[3] + 
                                                                        EPS*fem->psi[i][p]*F21_stiff2[q][3] + 
                                                                        EPS*fem->psi[i][q]*F21_stiff3[p][3] + 
                                                                         EPS*G_dual[1][0]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+4][q*this->dof_per_node+6] = 0.0;
                K[p*this->dof_per_node+4][q*this->dof_per_node+7] = 0.0;

                //********************** F_22 ***********************
                K[p*this->dof_per_node+5][q*this->dof_per_node] = 0.0;
                K[p*this->dof_per_node+5][q*this->dof_per_node+1] = -PEN*fem->psi[i][p]*(fem->dpsi[i][q][0]*G_dual[0][1] + 
                                                                                fem->dpsi[i][q][1]*G_dual[1][1]);
                K[p*this->dof_per_node+5][q*this->dof_per_node+2] = 0.0;
                K[p*this->dof_per_node+5][q*this->dof_per_node+3] = 0.0;
                K[p*this->dof_per_node+5][q*this->dof_per_node+4] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[0][1]) + 
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F22_stiff1[2] +
                                                                        EPS*fem->psi[i][p]*F22_stiff2[q][2] +
                                                                        EPS*fem->psi[i][q]*F22_stiff3[p][2] +  
                                                                         EPS*G_dual[0][1]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+5][q*this->dof_per_node+5] = PEN*fem->psi[i][p]*(fem->psi[i][q]*G_dual[1][1]) +
                                                                        EPS*fem->psi[i][p]*fem->psi[i][q]*F22_stiff1[3] + 
                                                                        EPS*fem->psi[i][p]*F22_stiff2[q][3] + 
                                                                        EPS*fem->psi[i][q]*F22_stiff3[p][3] + 
                                                                         EPS*G_dual[1][1]*
                                                                    (fem->dpsi[i][p][0]*fem->dpsi[i][q][0]*G_dual[0][0] + 
                                                                     fem->dpsi[i][p][0]*fem->dpsi[i][q][1]*G_dual[1][0] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][0]*G_dual[0][1] +
                                                                     fem->dpsi[i][p][1]*fem->dpsi[i][q][1]*G_dual[1][1]);

                K[p*this->dof_per_node+5][q*this->dof_per_node+6] = 0.0;
                K[p*this->dof_per_node+5][q*this->dof_per_node+7] = 0.0;

                //********************* beta_1 **********************
                K[p*this->dof_per_node+6][q*this->dof_per_node] = fem->psi[i][p]*(d_mu1_F[0]*fem->dpsi[i][q][0] + 
                                                                    d_mu1_F[1]*fem->dpsi[i][q][1]);

                K[p*this->dof_per_node+6][q*this->dof_per_node+1] = fem->psi[i][p]*(d_mu1_F[2]*fem->dpsi[i][q][0] + 
                                                                    d_mu1_F[3]*fem->dpsi[i][q][1]);

                K[p*this->dof_per_node+6][q*this->dof_per_node+2] = 0.0;
                K[p*this->dof_per_node+6][q*this->dof_per_node+3] = 0.0;
                K[p*this->dof_per_node+6][q*this->dof_per_node+4] = 0.0;
                K[p*this->dof_per_node+6][q*this->dof_per_node+5] = 0.0;
                K[p*this->dof_per_node+6][q*this->dof_per_node+6] = 0.0;
                K[p*this->dof_per_node+6][q*this->dof_per_node+7] = 0.0;

                //*********************** beta_2 **********************
                K[p*this->dof_per_node+7][q*this->dof_per_node] = fem->psi[i][p]*(d_mu2_F[0]*fem->dpsi[i][q][0] + 
                                                                    d_mu2_F[1]*fem->dpsi[i][q][1]);

                K[p*this->dof_per_node+7][q*this->dof_per_node+1] = fem->psi[i][p]*(d_mu2_F[2]*fem->dpsi[i][q][0] + 
                                                                    d_mu2_F[3]*fem->dpsi[i][q][1]);

                K[p*this->dof_per_node+7][q*this->dof_per_node+2] = 0.0;
                K[p*this->dof_per_node+7][q*this->dof_per_node+3] = 0.0;
                K[p*this->dof_per_node+7][q*this->dof_per_node+4] = 0.0;
                K[p*this->dof_per_node+7][q*this->dof_per_node+5] = 0.0;
                K[p*this->dof_per_node+7][q*this->dof_per_node+6] = 0.0;
                K[p*this->dof_per_node+7][q*this->dof_per_node+7] = 0.0;

            }
            
            // residual vector calculation
            // x_1
            fe[p*this->dof_per_node] = beta_1_gp*(d_mu1_F[0]*fem->dpsi[i][p][0] + d_mu1_F[1]*fem->dpsi[i][p][1]) + 
                                        beta_2_gp*(d_mu2_F[0]*fem->dpsi[i][p][0] + d_mu2_F[1]*fem->dpsi[i][p][1])

                                        - PEN*(F[0][0]*(G_dual[0][0]*fem->dpsi[i][p][0] + G_dual[0][1]*fem->dpsi[i][p][1])
                                            + F[0][1]*(G_dual[1][0]*fem->dpsi[i][p][0] + G_dual[1][1]*fem->dpsi[i][p][1]))

                                        + PEN*(dx1_zeta_1*fem->dpsi[i][p][0]*G_dual[0][0] + 
                                               dx1_zeta_1*fem->dpsi[i][p][1]*G_dual[0][1] +
                                               dx1_zeta_2*fem->dpsi[i][p][0]*G_dual[1][0] + 
                                               dx1_zeta_2*fem->dpsi[i][p][1]*G_dual[1][1]);
            // x_2
            fe[p*this->dof_per_node + 1] = beta_1_gp*(d_mu1_F[2]*fem->dpsi[i][p][0] + d_mu1_F[3]*fem->dpsi[i][p][1]) + 
                                            beta_2_gp*(d_mu2_F[2]*fem->dpsi[i][p][0] + d_mu2_F[3]*fem->dpsi[i][p][1])

                                            - PEN*(F[1][0]*(G_dual[0][0]*fem->dpsi[i][p][0] + G_dual[0][1]*fem->dpsi[i][p][1])
                                                + F[1][1]*(G_dual[1][0]*fem->dpsi[i][p][0] + G_dual[1][1]*fem->dpsi[i][p][1]))

                                            + PEN*(dx2_zeta_1*fem->dpsi[i][p][0]*G_dual[0][0] + 
                                               dx2_zeta_1*fem->dpsi[i][p][1]*G_dual[0][1] +
                                               dx2_zeta_2*fem->dpsi[i][p][0]*G_dual[1][0] + 
                                               dx2_zeta_2*fem->dpsi[i][p][1]*G_dual[1][1]);
            // F_11
            fe[p*this->dof_per_node + 2] = PEN*(fem->psi[i][p]*(F[0][0]*G_dual[0][0] + F[0][1]*G_dual[1][0]))
                                            - PEN*(fem->psi[i][p]*(dx1_zeta_1*G_dual[0][0] + dx1_zeta_2*G_dual[1][0]))
                                            + EPS*(F11_res1[p] + F11_res2[p] + F11_res3[p] + F11_res4[p]);
            // F_12
            fe[p*this->dof_per_node + 3] = PEN*(fem->psi[i][p]*(F[0][0]*G_dual[0][1] + F[0][1]*G_dual[1][1]))
                                             - PEN*(fem->psi[i][p]*(dx1_zeta_1*G_dual[0][1] + dx1_zeta_2*G_dual[1][1]))
                                             + EPS*(F12_res1[p] + F12_res2[p] + F12_res3[p] + F12_res4[p]);
            // F_21
            fe[p*this->dof_per_node + 4] = PEN*(fem->psi[i][p]*(F[1][0]*G_dual[0][0] + F[1][1]*G_dual[1][0]))
                                            - PEN*(fem->psi[i][p]*(dx2_zeta_1*G_dual[0][0] + dx2_zeta_2*G_dual[1][0]))
                                            + EPS*(F21_res1[p] + F21_res2[p] + F21_res3[p] + F21_res4[p]);
            // F_22
            fe[p*this->dof_per_node + 5] = PEN*(fem->psi[i][p]*(F[1][0]*G_dual[0][1] + F[1][1]*G_dual[1][1]))
                                            - PEN*(fem->psi[i][p]*(dx2_zeta_1*G_dual[0][1] + dx2_zeta_2*G_dual[1][1]))
                                            + EPS*(F22_res1[p] + F22_res2[p] + F22_res3[p] + F22_res4[p]);
            // beta_1
            fe[p*this->dof_per_node + 6] = fem->psi[i][p]*(mu1 - this->lambda_1*this->lambda_1);
            // beta_2
            fe[p*this->dof_per_node + 7] = fem->psi[i][p]*(mu2 - this->lambda_2*this->lambda_2);
        
        }


        for (int p=0;p<this->dof_per_el;p++){
			for (int q=0;q<this->dof_per_el;q++){

				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[i]*fem->wt[i];
			}

			rhs_el[p] = rhs_el[p] + (-1.0)*fe[p]*detJ_e[i]*fem->wt[i];
		}

    }// loop over gpts ends


    return 0;
    
}

PetscErrorCode Primal_Solve::deformation_gradient_l2_projection()
{

	PetscErrorCode ierr;

	// scalars
	PetscReal rnorm=0.0;
  	PetscInt  its=0;

	int dof_Fdef = fem->node_per_el*dual_solve->F_DOF;

	int indc[dof_Fdef];

	double **ke;

	ke = new double* [dof_Fdef];
	for (int i=0;i<dof_Fdef;i++)
		ke[i]= new double [dof_Fdef];

	double ke_vec[dof_Fdef*dof_Fdef];

	double fe[dof_Fdef];

	MyVec <double> detJ_e(fem->ngp);

    MyMat <double> x_def_el(fem->node_per_el,MyVec <double> (dual_solve->F_DOF));

	ierr = VecSet(F_fdef,0.0);CHKERRQ(ierr);

	int node_el;

	for (int ie=Istart;ie<Iend;++ie)
	{	

		// initialising to zero element stiff matrices and force vectors
		for (int k1=0;k1<dof_Fdef;k1++){
			for (int k2=0;k2<dof_Fdef;k2++){

				ke[k1][k2] = 0.0;
			}
			fe[k1] = 0.0;
		}

		// collecting jacobian for the current element
		for (int i=0;i<fem->ngp;i++){			
			detJ_e[i] = fem->detJ_G[ie][i];	
		}

        // collecting values at nodes of each element
		for (int i=0;i<fem->node_per_el;i++){

            node_el = grid->el_conn[ie][i];

			for (int k=0;k<dual_solve->F_DOF;k++){

				x_def_el[i][k] = x_def[node_el][k];
			}
		}

		ierr = this->elem_deformation_gradient_l2_projection_stiff_force(x_def_el,detJ_e,ke,fe); CHKERRQ(ierr);


		for (int k1=0;k1<fem->node_per_el;k1++){
			for (int d=0;d<dual_solve->F_DOF;d++){
				indc[dual_solve->F_DOF*k1+d] = dual_solve->F_DOF*grid->el_conn[ie][k1] + d;
			}
		}

		// adding values to global matrix and vector
		for (int i=0;i<dof_Fdef;i++){
			for (int j=0;j<dof_Fdef;j++){

				ke_vec[i*dof_Fdef+j] = ke[i][j];
			}
		}

		// adding values to global matrix and vector
		if (this->count == 0){
			ierr = MatSetValues(K_fdef,dof_Fdef,indc,dof_Fdef,indc,ke_vec,ADD_VALUES); CHKERRQ(ierr);
		}
		
		ierr = VecSetValues(F_fdef,dof_Fdef,indc,fe,ADD_VALUES); CHKERRQ(ierr);

	} // loop over the elements

	// final assembly of global vector
	ierr = VecAssemblyBegin(F_fdef); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_fdef); CHKERRQ(ierr);

	// final assembly of global matrix
	if (this->count == 0){
		ierr = MatAssemblyBegin(K_fdef,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(K_fdef,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}

	// creating the solver context
	ierr = common_utilities->ksp_mumps_solver_petsc(ksp_fdef,K_fdef,pc_fdef,FF_fdef); CHKERRQ(ierr);
	
	//Solve the linear system
	ierr = KSPSolve(ksp_fdef,F_fdef,fdef_vec);CHKERRQ(ierr);

	ierr = KSPGetResidualNorm(ksp_fdef,&rnorm); CHKERRQ(ierr);
	ierr = KSPGetTotalIterations(ksp_fdef,&its); CHKERRQ(ierr);

	if (::rank==0){
		std::cout << "deformation gradient solve with norm as " << rnorm << " and iterations are " << its << std::endl;
	}

	// context for scattering data to all processors
	if (this->count == 0){
		ierr = VecScatterCreateToAll(fdef_vec,&ctx_fdef,&fdef_vec_SEQ); CHKERRQ(ierr);
	}

	ierr = VecScatterBegin(ctx_fdef,fdef_vec,fdef_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx_fdef,fdef_vec,fdef_vec_SEQ,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = VecGetValues(fdef_vec_SEQ,tdof_fdef,ind_set_fdef,f_def_sol); CHKERRQ(ierr);


    // update the values of F_def on all nodes
	for (int i=0;i<grid->ndof;i++){
		for (int d=0;d<dual_solve->F_DOF;d++){

			F_def[i][d] = f_def_sol[dual_solve->F_DOF*i+d];
		}
	}


	// delete dynamically alotted matrix
	for (int i = 0; i<dof_Fdef; i++)
		delete[] ke[i];

	delete[] ke;


	return 0;

}

PetscErrorCode Primal_Solve::elem_deformation_gradient_l2_projection_stiff_force(MyMat <double> xdef_el,
                                                                            MyVec <double> detJ_e, double **ke, double *fe)
{   

    MyMat <double> F_gp(fem->ngp,MyVec <double> (dual_solve->F_DOF));

    int dof_fdef = fem->node_per_el*dual_solve->F_DOF;

	double K[dof_fdef][dof_fdef];
	double F[dof_fdef];

	// getting F_gp from x_def_el at all gps
	for(int i=0;i<fem->ngp; i++)
    {

        F_gp[i][0] = 0.0; F_gp[i][1] = 0.0; F_gp[i][2] = 0.0; F_gp[i][3] = 0.0;

        for (int p=0;p<fem->node_per_el;p++){
            
            F_gp[i][0] = F_gp[i][0] + fem->dpsi[i][p][0]*xdef_el[p][0];
            F_gp[i][1] = F_gp[i][1] + fem->dpsi[i][p][1]*xdef_el[p][0];
            F_gp[i][2] = F_gp[i][2] + fem->dpsi[i][p][0]*xdef_el[p][1];
            F_gp[i][3] = F_gp[i][3] + fem->dpsi[i][p][1]*xdef_el[p][1];
        }
    }

	for (int gp=0;gp<fem->ngp;gp++)
	{

        // initializing to zero
		for (int p=0;p<dof_fdef;p++){
			for (int q=0;q<dof_fdef;q++){
				K[p][q] = 0.0;
			}
			F[p] = 0.0;
		}


		for (int p=0;p<fem->node_per_el;p++){
			for (int q=0;q<fem->node_per_el;q++){

				for (int d=0;d<dual_solve->F_DOF;d++)
				{
                    // all other terms are zero
					K[dual_solve->F_DOF*p+d][dual_solve->F_DOF*q+d] = fem->psi[gp][p]*fem->psi[gp][q];	
				}// stiffness matrix
				
			}

			for (int d=0;d<dual_solve->F_DOF;d++){
				F[dual_solve->F_DOF*p + d] = F_gp[gp][d]*fem->psi[gp][p];
			}// force vector
			
		}

        // summation at gauss point
		for(int p=0;p<dof_fdef;p++)
        {
			for(int q=0;q<dof_fdef;q++)
            {
				ke[p][q] = ke[p][q] + K[p][q]*detJ_e[gp]*fem->wt[gp];
			}

			fe[p] = fe[p] + F[p]*detJ_e[gp]*fem->wt[gp];
		}

	}// loop over the gauss pts

	return 0;
}
