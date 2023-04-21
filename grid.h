#ifndef GRID_H
#define GRID_H

class Grid
{
private:

public:
	Grid(int nx, int ny, double l1, double l2);

	~Grid();

	int nx,ny;
	double l1,l2;

	int ndof; // no of nodes
	int nel; // no of elements in 2d
	int nbel; // no of boundary nodes and elements
	int tdof; // total dof for dual functional
	double h_min; // min h in the mesh used for time step
	int dim;
	int A_DOF;
	int BETA_DOF;
	int X_DOF;
	int dof_per_node;

	// co-ordinates
	MyMat<double> X_ref;

	// el-connectivity (1 dof)
	MyMat<int> el_conn;

	// el-connectivity (6 or 8 dof) for dual
	MyMat<int> dof_el_conn;

	// boundary el-conn (1 dof) 
	MyMat<int> b_el_conn;

	// boundary el-conn (6 or 8 dof) for dual
	MyMat<int> dof_b_el_conn;

	// boundary ids
	MyVec<int> boundary_ids;

	void set_mesh_parameters();

	void mesh_info_file();
	
	void initialize_fields();

	PetscErrorCode mesh_data_read();

	PetscErrorCode mesh_rectangular();

	PetscErrorCode mesh_trelis();

	PetscErrorCode get_h_min();

};


#endif