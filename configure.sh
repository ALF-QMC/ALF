export DIR=`pwd`
export f90=mpiifort
#export f90="gfortran"
export FL="-g -c  -w  -O3 -xHost -unroll -finline-functions -heap-arrays 1024"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-mkl=sequential"
# uncomment next lines for hybrid of MPI and OpenMP
export LIB_BLAS_LAPACK="-mkl=parallel -qopenmp"
export FL="${FL} -qopenmp -parallel"
