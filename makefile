-include ../petscdir.mk
CFLAGS           =
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/sys/tutorials/
DIRS             = network
EXAMPLESC        = ex12.c solver.c main_code.c
EXAMPLESF        = 
MANSEC           = Sys
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
	
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
#include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

ex12: ex12.o  chkopts
	-${CLINKER} -o ex12 ex12.o  ${PETSC_KSP_LIB}
	${RM} ex12.o

solver: solver.o  chkopts 
	-${CLINKER} -o solver solver.o  ${PETSC_KSP_LIB}
	${RM} solver.o

lcg_problem: lcg_problem.o 
	-${CLINKER} -o lcg_problem lcg_problem.o ${PETSC_KSP_LIB}
	${RM} lcg_problem.o
	
main: main.o chkopts 
	-${CLINKER} -o main main.o  ${PETSC_KSP_LIB}
	${RM} main.o

#clean:
#	rm *.o rupture_layer
