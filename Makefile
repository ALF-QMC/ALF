# -DMPI selects MPI.
# -DSTAB1  Alternative stabilization, using the singular value decomposition.
# -DSTAB2  Alternative stabilization, lapack QR with  manual pivoting. Packed form of QR factorization is not used.
# (Noflag) Default  stabilization, using lapack QR with pivoting. Packed form of QR factorization  is used. 
# -DQRREF  Enables reference lapack implementation of QR decomposition.
# -DTEMPERING  Complies program for parallel tempering. Requires MPI
# Recommendation:  just use the -DMPI flag if you want to run in parallel or leave it empy for serial jobs.  
# The default stabilization, no flag, is generically the best. 
PROGRAMCONFIGURATION = -DMPI 
PROGRAMCONFIGURATION = 
PROGRAMCONFIGURATION = -DMPI  -DTEMPERING
f90 = gfortran
f90 = $(mpif90)
export f90
F90OPTFLAGS = -O3 -Wconversion  -fcheck=all
F90OPTFLAGS = -O3
export F90OPTFLAGS
F90USEFULFLAGS = -cpp -std=f2003
F90USEFULFLAGS = -cpp
export F90USEFULFLAGS
FL = -c ${F90OPTFLAGS} ${PROGRAMCONFIGURATION}
export FL
DIR = ${CURDIR}
export DIR
Libs = ${DIR}/Libraries/
export Libs
LIB_BLAS_LAPACK = -llapack -lblas
export LIB_BLAS_LAPACK

all: lib ana program  Hub_Ising SPT Hub Hub_Can Kondo_Honey

lib:
	cd Libraries && $(MAKE)
ana:
	cd Analysis && $(MAKE)
program:
	cd Prog && $(MAKE) Examples
Hub_Ising:
	cd Prog && $(MAKE) Hub_Ising
SPT:
	cd Prog && $(MAKE) SPT
Hub:
	cd Prog && $(MAKE) Hub
Hub_Can:
	cd Prog && $(MAKE) Hub_Can
Kondo_Honey:
	cd Prog && $(MAKE) Kondo_Honey


clean: cleanall
cleanall: cleanprog cleanlib cleanana  
cleanprog:
	cd Prog && $(MAKE) clean 
cleanlib:
	cd Libraries && $(MAKE) clean
cleanana:
	cd Analysis && $(MAKE) clean
help:
	@echo "The following are some of the valid targets of this Makefile"
	@echo "all, program, lib, ana, clean, cleanall, cleanprog, cleanlib, cleanana"
	@echo "Hub_Ising SPT Hub Hub_Can Kondo_Honey"
