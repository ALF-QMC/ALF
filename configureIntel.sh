export DIR=`pwd`
export f90=mpiifort
#export f90="gfortran"

export FL="-c -O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024"

#uncomment the next line if you want to debug/profile your code
#export FL="${FL} -g -traceback"

export Libs=${DIR}"/Libraries/"
#export LIB_BLAS_LAPACK="-lblas -llapack"
export LIB_BLAS_LAPACK="-mkl=sequential"
