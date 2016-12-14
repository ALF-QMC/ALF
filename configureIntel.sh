export DIR=`pwd`
export f90=ifort
#export f90="gfortran"

export FL="-c -O3 -fp-model fast=2 -axCORE-AVX2 -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin -DnMPI -DZLIB"
# -parallel -qopenmp"

#uncomment the next line if you want to debug/profile your code
export FL="${FL} -g -traceback"

export Libs=${DIR}"/Libraries/"
#export LIB_BLAS_LAPACK="-lblas -llapack"
export LIB_BLAS_LAPACK="-mkl -lz"
#export LIB_BLAS_LAPACK="-mkl=parallel"
