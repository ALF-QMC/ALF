export DIR=`pwd`
export f90=mpiifort
#export f90="gfortran"

#LRZ enviroment
export f90=mpif90
module unload mpi.ibm
module load mpi.intel
module unload mkl
module load mkl/2017_s

export FL="-c  -O3  -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024"
#uncomment the next line if you want to debug/profile your code
export FL="${FL} -g -traceback"

export Libs=${DIR}"/Libraries/"
#export LIB_BLAS_LAPACK="-lblas -llapack"
#export LIB_BLAS_LAPACK="-mkl=sequential"

#uncomment the following line on LRZ
export LIB_BLAS_LAPACK=$MKL_LIB

# uncomment next lines for hybrid of MPI and OpenMP
#export LIB_BLAS_LAPACK="-mkl=parallel -qopenmp"
#export FL="${FL} -qopenmp -parallel"
