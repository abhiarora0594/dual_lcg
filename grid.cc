#include "grid.h"

Grid::Grid(int nx, int ny, double l1, double l2)
{	

	this->nx = nx;
	this->ny = ny;
	this->l1 = l1;
	this->l2 = l2;

	if (problem_type == inverse){
		this->dim = 2;
	}
	else if (problem_type == forward){
		this->dim = 3;
	}

	this->A_DOF = 2*(this->dim);
	this->BETA_DOF = 2;
	this->X_DOF = 3;
	this->dof_per_node = this->BETA_DOF + this->A_DOF;

	this->h_min = this->l1;

	// setting mesh variables
	if (mesh_input == cylindrical || mesh_input == plane){
		this->set_mesh_parameters();
	}
	else if (mesh_input == trelis_mesh || mesh_input == spherical){
		this->mesh_info_file();
	}

	// allocating space for fields
	this->initialize_fields();

}

Grid::~Grid()
{

}

void Grid::set_mesh_parameters()
{
	ndof = (nx+1)*(ny+1);
	nel = nx*ny;
	nbel = 2*nx + 2*ny;
	tdof = (BETA_DOF+A_DOF)*ndof;

}

void Grid::mesh_info_file()
{
	std::ifstream myfile;

	myfile.open("./input/mesh_info.in",std::ios::in);

	myfile >> ndof >> nel >> nbel;

	tdof = (BETA_DOF+A_DOF)*ndof; // (2+4 or 2+6) dof at each node

	myfile.close();

}

void Grid::initialize_fields()
{

	// co-ordinates
	X_ref = MyMat<double> (ndof,MyVec<double>(X_DOF));

	// el-connectivity (1 dof)
	el_conn = MyMat<int> (ndof,MyVec<int>(4));

	// el-connectivity (6 or 8 dof) for dual
	dof_el_conn = MyMat<int> (ndof,MyVec<int>(4*(BETA_DOF+A_DOF)));

	// boundary el-conn (1 dof) 
	b_el_conn = MyMat<int> (nbel,MyVec<int>(2));

	// boundary el-conn (6 or 8 dof) for dual
	dof_b_el_conn = MyMat<int> (nbel,MyVec<int>(2*(BETA_DOF+A_DOF)));

	// boundary ids
	boundary_ids = MyVec<int> (nbel);

}


PetscErrorCode Grid::mesh_data_read()
{

	int count;

	PetscErrorCode ierr;

	std::cout << "ndof is " << ndof << std::endl;

	std::cout << "size of X_ref is " << X_ref.size() << " and " << X_ref[0].size() << std::endl;


	if (mesh_input == cylindrical || mesh_input == plane){
		ierr = this->mesh_rectangular(); CHKERRQ(ierr);
	}
	else if(mesh_input == trelis_mesh || mesh_input == spherical){
		ierr = this->mesh_trelis(); CHKERRQ(ierr);
	}

	for (int i=0;i<nel;i++){

		count = 0;

		// element connectivity for (6 or 8 dof per node) field
		for (int k=0;k<(BETA_DOF+A_DOF);k++){
			dof_el_conn[i][count] = (BETA_DOF+A_DOF)*el_conn[i][0]+k;
			count = count + 1; 
		}

		for (int k=0;k<(BETA_DOF+A_DOF);k++){
			dof_el_conn[i][count] = (BETA_DOF+A_DOF)*el_conn[i][1]+k; 
			count = count + 1;
		}

		for (int k=0;k<(BETA_DOF+A_DOF);k++){
			dof_el_conn[i][count] = (BETA_DOF+A_DOF)*el_conn[i][2]+k; 
			count = count + 1;
		}

		for (int k=0;k<(BETA_DOF+A_DOF);k++){
			dof_el_conn[i][count] = (BETA_DOF+A_DOF)*el_conn[i][3]+k; 
			count = count + 1;
		}

	}


	for (int i=0;i<nbel;i++){

		count=0;

		for (int k=0;k<(BETA_DOF+A_DOF);k++){
			dof_b_el_conn[i][count] = (BETA_DOF+A_DOF)*b_el_conn[i][0]+k;
			count = count + 1;
		}

		for (int k=0;k<(BETA_DOF+A_DOF);k++){
			dof_b_el_conn[i][count] = (BETA_DOF+A_DOF)*b_el_conn[i][1]+k;
			count = count + 1;
		}
	}

	return 0;

}

PetscErrorCode Grid::mesh_rectangular()
{

	double r=l1;
	double z=l2;
	int tmp;
	double theta;

	if (mesh_input == cylindrical)
	{
		for (int i=0;i<(nx+1);++i){
			for (int j=0;j<(ny+1);++j){

				tmp = i + j*(nx+1); // node no
				theta = PI*(1.0 - ((double)i/(double)nx));
				X_ref[tmp][0] = r*std::cos(theta); // x-coordinate
				X_ref[tmp][2] = r*std::sin(theta); // z-coordinate
				X_ref[tmp][1] = -z/2.0 + z*((double)j/(double)ny); // y-coordinate
			}
		}
	}
	else if (mesh_input == plane)
	{
		for (int i=0;i<(nx+1);++i){
			for (int j=0;j<(ny+1);++j){

				tmp = i + j*(nx+1); // node no
				X_ref[tmp][0] = -r/2.0 + r*((double)i/(double)nx); // x-coordinate
				X_ref[tmp][2] = 0.0; // z-coordinate
				X_ref[tmp][1] = -z/2.0 +  z*((double)j/(double)ny); // y-coordinate
			}
		}
	}

	for (int i=0;i<nx;++i){
		for (int j=0;j<ny;++j){

			tmp = i + j*nx; // element no

			// element connectivity for (1 dof per node) field
			el_conn[tmp][0] = (i+j*(nx+1));
			el_conn[tmp][1] = (i+1+j*(nx+1)); 
			el_conn[tmp][2] = (i+1+(j+1)*(nx+1)); 
			el_conn[tmp][3] = (i+(j+1)*(nx+1)); 

		}
	}

	int count=0;

	// bottom boundary
	for (int i=0;i<nx;i++){

		b_el_conn[count][0]	= i; 
		b_el_conn[count][1] = i+1; 
		boundary_ids[count] = 1; // bottom

		count = count+1;
	}

	// right boundary
	for (int i=0;i<ny;i++){

		b_el_conn[count][0]	= (i+1)*(nx+1)-1;
		b_el_conn[count][1] = (i+2)*(nx+1)-1;	
		boundary_ids[count] = 2; // right

		count = count+1;
	}

	// top boundary
	for (int i=0;i<nx;i++){

		b_el_conn[count][0]	= (nx+1)*(ny+1)-i-1; 
		b_el_conn[count][1] = (nx+1)*(ny+1)-i-2;
		boundary_ids[count] = 3; // top

		count = count+1;
	}

	// left boundary
	for (int i=0;i<ny;i++){

		b_el_conn[count][0]	= (nx+1)*(ny-i); 
		b_el_conn[count][1] = (nx+1)*(ny-i-1);
		boundary_ids[count] = 4; // left 

		count = count+1;
	}

	return 0;

}

PetscErrorCode Grid::mesh_trelis()
{

	double val;

	std::ifstream myfile1,myfile2,myfile3;

	myfile1.open("./input/coordinates.in",std::ios::in);

	for (int i=0;i<ndof;++i){

		myfile1 >> X_ref[i][0] >> X_ref[i][1] >> X_ref[i][2];
	}

	myfile1.close();

	if (mesh_input == spherical){

		double r = 1.0;
		double r_0 = -5.0; //-0.95
		double R = std::sqrt(r*r + r_0*r_0);
		double z_c = R + ((r_0-R)*2.0*R*(R-r_0)/((R-r_0)*(R-r_0) + r*r));

		for (int i=0;i<ndof;++i){

			val = (R-r_0)*(R-r_0) + X_ref[i][0]*X_ref[i][0] + X_ref[i][1]*X_ref[i][1];
			X_ref[i][0] = X_ref[i][0];//2.0*R*(R-r_0)*X_ref[i][0]/val;
			X_ref[i][1] = X_ref[i][1];//2.0*R*(R-r_0)*X_ref[i][1]/val;
			X_ref[i][2] = 0.0;// R + ((r_0-R)*2.0*R*(R-r_0)/val) - z_c;
		}
	}

	myfile2.open("./input/element_connectivity.in",std::ios::in);

	for (int i=0;i<nel;++i){

		myfile2 >> el_conn[i][0] >> el_conn[i][1] >> el_conn[i][2] >> el_conn[i][3];
	}

	myfile2.close();

	myfile3.open("./input/boundary_elems.in",std::ios::in);

	for (int i=0;i<nbel;++i){
		myfile3 >> b_el_conn[i][0] >> b_el_conn[i][1];
	}

	myfile3.close();

	for (int i=0;i<nbel;++i){
		boundary_ids[i] = 0; // for any random shape
	}

	return 0;

}


PetscErrorCode Grid::get_h_min()
{

	double p1[3], p2[3], p3[3];

	double val1, val2, val3;

	for (int i=0;i<nel;i++){

		p1[0] = X_ref[el_conn[i][0]][0];
		p1[1] = X_ref[el_conn[i][0]][1];
		p1[2] = X_ref[el_conn[i][0]][2];

		p2[0] = X_ref[el_conn[i][1]][0];
		p2[1] = X_ref[el_conn[i][1]][1];
		p2[2] = X_ref[el_conn[i][1]][2];

		p3[0] = X_ref[el_conn[i][3]][0];
		p3[1] = X_ref[el_conn[i][3]][1];
		p3[2] = X_ref[el_conn[i][3]][2];

		val1 = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]));
		val2 = std::sqrt((p1[0]-p3[0])*(p1[0]-p3[0]) + (p1[1]-p3[1])*(p1[1]-p3[1]) + (p1[2]-p3[2])*(p1[2]-p3[2]));

		val3 = std::min(val1,val2);

		this->h_min = std::min(this->h_min,val3);
	}

	std::cout << "h_min is " << this->h_min << std::endl;

	return 0;
}