export DIR=`pwd`

#LRZ enviroment
module switch mpi.ibm mpi.intel
module switch intel intel/17.0
module switch mkl mkl/2017
export f90=mpif90

export FL="-c -O3 -fp-model fast=2 -xCORE-AVX2 -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin -DMPI -DZLIB"
# -parallel -qopenmp"
#uncomment the next line if you want to debug/profile your code
#export FL="${FL} -g -traceback"

export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK=$MKL_LIB" -lz"
