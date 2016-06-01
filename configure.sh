export DIR=`pwd`
export f90=mpiifort
#export f90="gfortran"
export FL="-c  -w  -O3 -heap-arrays 1024"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-mkl=sequential"
