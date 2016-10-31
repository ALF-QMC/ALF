export DIR=`pwd`
export f90=ifort
#export f90="gfortran"

export FL="-c -O3 -fp-model fast=2 -march=corei7-avx -unroll -finline-functions -ipo -ip -heap-arrays 1024 -DnMPI -parallel -qopenmp"

#uncomment the next line if you want to debug/profile your code
export FL="${FL} -g"

export Libs=${DIR}"/Libraries/"
#export LIB_BLAS_LAPACK="-lblas -llapack"
#export LIB_BLAS_LAPACK="-mkl=sequential"
export LIB_BLAS_LAPACK="-mkl=parallel"
