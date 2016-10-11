export DIR=`pwd`

#LRZ enviroment
export f90=mpif90
module unload mpi.ibm
module load mpi.intel
module unload mkl
module load mkl/11.3_s

export FL="-c -O3 -std03 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024"
#uncomment the next line if you want to debug/profile your code
export FL="${FL} -g -traceback"

export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK=$MKL_LIB
