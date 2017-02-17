export DIR=`pwd`
export f90=ifort
#export f90="gfortran"

export FL="-c -O3 -fp-model fast=2 -axCORE-AVX2 -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin -DnMPI"
#uncomment the next line if you want compressed tau-resolved data files
#export FL="${FL} -DZLIB"

#export FL="-c -O0 -heap-arrays 1024 -no-wrap-margin"
# -parallel -qopenmp"

#uncomment the next line if you want to debug/profile your code
export FL="${FL} -g -traceback"

export Libs=${DIR}"/Libraries/"
#export LIB_BLAS_LAPACK="-lblas -llapack"
export LIB_BLAS_LAPACK="-mkl -lz"
#export LIB_BLAS_LAPACK="-mkl=parallel"
