#ifndef FEM_H
#define FEM_H

class Grid;

class Fem
{
private:
	void gauss_pts();

	void local_sf_der_linear();

public:
	Fem(Grid *grid, int ngp = 4, int node_per_el = 4);

	~Fem();

	int ngp;

	int node_per_el;

	int dof_per_el;

	Grid *grid;
	
	MyMat <double> xi;

	MyVec <double> wt;

	MyMat <double> psi;

	MyMat <MyVec<double>> dpsi;

	MyVec <double> psi_c;

	MyMat <double> dpsi_c;

	MyMat <MyVec<double>> E_real;

	MyMat <MyVec<double>> E_dual;

	MyMat <double> E_n_gp;

	MyMat <double> detJ_G;

	PetscErrorCode jacobian_convected_basis();

	template <class T>
	PetscErrorCode dot_product(MyVec<T> v1, MyVec<T> v2, T *tmp);

	template<class T>
	PetscErrorCode my_mat_inverse(MyMat<T> g, MyMat<T> & g_inv);

	template<class T>
	PetscErrorCode cross_product(MyVec<T> v1, MyVec<T>v2, MyVec<T>& v3);

	template<class T>
	PetscErrorCode jacobian(MyVec<T> v1, MyVec<T> v2, T *tmp_sc);
	
};


#endif