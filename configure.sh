export DIR=`pwd`
export f90=mpiifort
#export f90="gfortran"

export FL="-c  -O3  -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024"
#uncomment the next line if you want to debug/profile your code
export FL="${FL} -g"

export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-lblas -llapack"
export LIB_BLAS_LAPACK="-mkl=sequential"
# uncomment next lines for hybrid of MPI and OpenMP
#export LIB_BLAS_LAPACK="-mkl=parallel -qopenmp"
#export FL="${FL} -qopenmp -parallel"
