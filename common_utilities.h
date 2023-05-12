#ifndef COMMON_UTILITIES_H
#define COMMON_UTILITIES_H

#include "dual_solve.h"

class Dual_Solve;
class Fem;
class Grid;

class Common_Utilities
{
private:

    PetscInt Istart,Iend;

public:
	Common_Utilities(Fem *fem, Grid *grid, Dual_Solve *dual_solve);

	~Common_Utilities();

    Fem *fem;
    Grid *grid;
    Dual_Solve *dual_solve;

    // methods
	PetscErrorCode vtk_write();

	PetscErrorCode restart_write();

	PetscErrorCode mysolve_local_linear_system(int , MyMat <double> ,
												MyVec <double> , MyVec <double> &);

	PetscErrorCode petsc_local_linear_system(int , MyMat <double> ,
												MyVec <double> , MyVec <double> &);

    PetscErrorCode mat_create_petsc(Mat &, int, int, int);

	PetscErrorCode vec_create_petsc(Vec &, int);

	PetscErrorCode ksp_mumps_solver_petsc(KSP &, Mat &, PC &, Mat &);

    template<class T>
	inline PetscErrorCode get_max_val(MyVec<T> v, T *val);

};


#endif