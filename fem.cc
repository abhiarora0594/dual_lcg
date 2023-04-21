#include "fem.h"


Fem::Fem(Grid *grid, int ngp, int node_per_el)
{
	
	this->grid = grid;
	this->ngp = ngp;
	this->node_per_el = node_per_el;
	this->dof_per_el = this->node_per_el*grid->dof_per_node;

	xi = MyMat <double> (ngp,MyVec<double>(2));

	wt = MyVec <double> (ngp);

	psi = MyMat <double> (ngp,MyVec<double>(4));

	dpsi = MyMat <MyVec<double>> (ngp,MyMat<double>(4,MyVec<double>(2)));

	psi_c = MyVec <double> (ngp);

	dpsi_c = MyMat <double> (ngp,MyVec<double>(2));

	E_real = MyMat <MyVec<double>> (grid->nel,MyMat<double>(ngp,MyVec<double>(2*grid->X_DOF)));

	E_dual = MyMat <MyVec<double>> (grid->nel,MyMat<double>(ngp,MyVec<double>(2*grid->X_DOF)));

	E_n_gp = MyMat <double> (grid->nel,MyVec<double>(ngp*grid->X_DOF));

	detJ_G = MyMat <double> (grid->nel,MyVec<double>(ngp));

	this->gauss_pts();

	this->local_sf_der_linear();
}

Fem::~Fem()
{

}

void Fem::gauss_pts()
{
	if (ngp==4){

		// 1st index is for gauss pt, 2nd is for (zeta_1 and zeta_2) co-ordinate
		xi[0][0] = -0.577350269189626; xi[0][1] = -0.577350269189626;
		xi[1][0] = 0.577350269189626; xi[1][1] = -0.577350269189626;
		xi[2][0] = 0.577350269189626; xi[2][1] = 0.577350269189626;
		xi[3][0] = -0.577350269189626; xi[3][1] = 0.577350269189626;

		for (int i=0;i<ngp;i++){
			wt[i] = 1.0;
		}
	}
}

void Fem::local_sf_der_linear()
{

	// psi (1st index is for gauss point, 2nd for shape fn)
	// dpsi (1st index is for gauss point, 2nd for shape fn, 3rd for zeta_1 or zeta_2 derivative)

	double sqrt_three = 0.577350269189626;

	for (int i=0;i<ngp;i++){

		psi[i][0] = 0.25*(1.0-xi[i][0])*(1.0-xi[i][1]);
		psi[i][1] = 0.25*(1.0+xi[i][0])*(1.0-xi[i][1]);
		psi[i][2] = 0.25*(1.0+xi[i][0])*(1.0+xi[i][1]);
		psi[i][3] = 0.25*(1.0-xi[i][0])*(1.0+xi[i][1]);

	}


	for (int i=0;i<ngp;i++){

		dpsi[i][0][0]= -0.25*(1.0-xi[i][1]);
		dpsi[i][1][0]= 0.25*(1.0-xi[i][1]);
		dpsi[i][2][0]= 0.25*(1.0+xi[i][1]);
		dpsi[i][3][0]= -0.25*(1.0+xi[i][1]);


		dpsi[i][0][1]= -0.25*(1.0-xi[i][0]);
		dpsi[i][1][1]= -0.25*(1.0+xi[i][0]);
		dpsi[i][2][1]= 0.25*(1.0+xi[i][0]);
		dpsi[i][3][1]= 0.25*(1.0-xi[i][0]);
	}



	// shape funcs at the center
	for (int i=0;i<ngp;i++){

		psi_c[0] = 0.25;
		psi_c[1] = 0.25;
		psi_c[2] = 0.25;
		psi_c[3] = 0.25;

	}

	// derivatives of shape funcs at the center
	for (int i=0;i<ngp;i++){

		dpsi_c[0][0]= -0.25*sqrt_three;
		dpsi_c[1][0]= 0.25*sqrt_three;
		dpsi_c[2][0]= 0.25*sqrt_three;
		dpsi_c[3][0]= -0.25*sqrt_three;


		dpsi_c[0][1]= -0.25*sqrt_three;
		dpsi_c[1][1]= -0.25*sqrt_three;
		dpsi_c[2][1]= 0.25*sqrt_three;
		dpsi_c[3][1]= 0.25*sqrt_three;
	}


	// for (int k1=0;k1<4;k1++){
	// 	for (int k2=0;k2<2;k2++){
	// 			std::cout << "dpsi is " << dpsi[0][k1][k2] << std::endl;			
	// 	}
	// }
		
}

PetscErrorCode Fem::jacobian_convected_basis()
{

	PetscErrorCode ierr;

	MyVec<double> tmp_1 = MyVec<double>(3);
	MyVec<double> tmp_2 = MyVec<double>(3);
	MyVec<double> tmp_3 = MyVec<double>(3);

	MyMat<double> G_ij = MyMat<double> (2,MyVec<double>(2));
	MyMat<double> G_inv = MyMat<double> (2,MyVec<double>(2));

	double x_el[4],y_el[4],z_el[4];
	double tmp_sc;

	for (int ie=0;ie<grid->nel;ie++)
	{

		for (int k=0;k<node_per_el;k++){
			x_el[k] = grid->X_ref[grid->el_conn[ie][k]][0];
			y_el[k] = grid->X_ref[grid->el_conn[ie][k]][1];
			z_el[k] = grid->X_ref[grid->el_conn[ie][k]][2];
		}

		for (int i=0;i<ngp;i++){

			// natural basis
			for (int alp=0;alp<2;++alp){

				E_real[ie][i][3*alp]= x_el[0]*dpsi[i][0][alp] + x_el[1]*dpsi[i][1][alp] + x_el[2]*dpsi[i][2][alp] + x_el[3]*dpsi[i][3][alp];
				E_real[ie][i][3*alp+1]= y_el[0]*dpsi[i][0][alp] + y_el[1]*dpsi[i][1][alp] + y_el[2]*dpsi[i][2][alp] + y_el[3]*dpsi[i][3][alp];
				E_real[ie][i][3*alp+2]= z_el[0]*dpsi[i][0][alp] + z_el[1]*dpsi[i][1][alp] + z_el[2]*dpsi[i][2][alp] + z_el[3]*dpsi[i][3][alp];
			}
			
			for (int alp=0;alp<2;++alp){
				for (int beta=0;beta<2;++beta){

					tmp_1[0] = E_real[ie][i][3*alp];
					tmp_1[1] = E_real[ie][i][3*alp+1];
					tmp_1[2] = E_real[ie][i][3*alp+2];

					tmp_2[0] = E_real[ie][i][3*beta];
					tmp_2[1] = E_real[ie][i][3*beta+1];
					tmp_2[2] = E_real[ie][i][3*beta+2];

					// convected basis
					ierr = this->dot_product(tmp_1,tmp_2,&tmp_sc); CHKERRQ(ierr);
					G_ij[alp][beta] = tmp_sc;

				}
			}

			// inverse of G_ij matrix
			ierr = this->my_mat_inverse(G_ij,G_inv); CHKERRQ(ierr);

			// dual basis 
			for (int alp=0;alp<2;++alp){

				E_dual[ie][i][3*alp]= G_inv[alp][0]*E_real[ie][i][0] + G_inv[alp][1]*E_real[ie][i][3];
				E_dual[ie][i][3*alp+1]= G_inv[alp][0]*E_real[ie][i][1] + G_inv[alp][1]*E_real[ie][i][4];
				E_dual[ie][i][3*alp+2]= G_inv[alp][0]*E_real[ie][i][2] + G_inv[alp][1]*E_real[ie][i][5];
			}

			// Jacobian of FE element
			tmp_1[0] = E_real[ie][i][0];
			tmp_1[1] = E_real[ie][i][1];
			tmp_1[2] = E_real[ie][i][2];

			tmp_2[0] = E_real[ie][i][3];
			tmp_2[1] = E_real[ie][i][4];
			tmp_2[2] = E_real[ie][i][5];

			ierr = this->jacobian(tmp_1,tmp_2,&tmp_sc); CHKERRQ(ierr);

			ierr = this->cross_product(tmp_1,tmp_2,tmp_3); CHKERRQ(ierr);

			// E_1 \cross E_2
			E_n_gp[ie][3*i] = tmp_3[0];
			E_n_gp[ie][3*i+1] = tmp_3[1];
			E_n_gp[ie][3*i+2] = tmp_3[2];

			detJ_G[ie][i] = tmp_sc; // jacobian


			// std::cout << "E_real for ie " << ie << " and gp " << i << " is " << E_real[ie][i][0] << " , "
			// 																	<< E_real[ie][i][1] << " , "
			// 																	<< E_real[ie][i][2] << " , "
			// 																	<< E_real[ie][i][3] << " , "
			// 																	<< E_real[ie][i][4] << " , "
			// 																	<< E_real[ie][i][5] << std::endl;

			// std::cout << "E_dual for ie " << ie << " and gp " << i << " is " << E_dual[ie][i][0] << " , "
			// 																	<< E_dual[ie][i][1] << " , "
			// 																	<< E_dual[ie][i][2] << " , "
			// 																	<< E_dual[ie][i][3] << " , "
			// 																	<< E_dual[ie][i][4] << " , "
			// 																	<< E_dual[ie][i][5] << std::endl;	

		}//loop over the gauss points

	}//loop over the elements

	std::cout << "size of E_real is " << E_real.size() << " and " <<
										 E_real[0].size() << " and " <<
										 E_real[0][0].size() << std::endl;

	return 0;
}

template<class T>
PetscErrorCode Fem::dot_product(MyVec<T> v1, MyVec<T> v2, T *tmp)
{

	T val = 0;
	// throw an exception if v1 and v2 size is not same

	for (unsigned int i=0;i<v1.size();i++){
		val = val + v1[i]*v2[i];
	}

	*tmp = val;

	return 0;
}


template<class T>
PetscErrorCode Fem::my_mat_inverse(MyMat<T> g, MyMat<T> & g_inv)
{	

	double detg;
	double detg_inv;

	if (g.size() == 2 && g[0].size()==2)
	{
		// determinant
		detg = g[1][1]*g[0][0] - g[0][1]*g[1][0];

		if (detg > 1e-12)
			detg_inv = 1.0/detg;
		else{
			// detg_inv = 0.0;
			std::cout << "determinant of the matrix is " << detg << std::endl;
			throw "determinant of matrix is 0 or negative, not invertible";
		}

			//detg_inv = 0.0;

		//std::cout << "Det is " << detg<< std::endl;

		g_inv[0][0] = detg_inv*g[1][1];
		g_inv[0][1] = -detg_inv*g[0][1];
		g_inv[1][0] = -detg_inv*g[1][0];
		g_inv[1][1] = detg_inv*g[0][0];	
	}
	// throw an exception if size of g is something else

	return 0;
}

template<class T>
PetscErrorCode Fem::cross_product(MyVec<T> v1, MyVec<T>v2, MyVec<T> & v3)
{

	if (v1.size()==3 && v2.size()==3){

		v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
		v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
		v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}
	else if (v1.size()==2 && v2.size()==2){
		
		v3[0] = 0.0;
		v3[1] = 0.0;
		v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}
	// throw an exception if v1 and v2 size is not same	

	return 0;
}

template<class T>
PetscErrorCode Fem::jacobian(MyVec<T> v1, MyVec<T> v2, T *tmp_sc)
{

	PetscErrorCode ierr;

	if (v1.size()==3 && v2.size()==3){
		
		MyVec<T> tmp = MyVec<T> (v1.size());

		ierr = this->cross_product(v1,v2,tmp); CHKERRQ(ierr);	

		*tmp_sc = std::sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
	}
	// throw an exception if v1 and v2 size is not same

	return 0;

}