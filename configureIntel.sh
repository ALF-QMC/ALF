export DIR=`pwd`
export f90=ifort
#export f90="gfortran"

PROGRAMMCONFIGURATION=""
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DSTAB1"
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DSTAB2"
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DQRREF"
# uncomment the next line if you want an MPI parallel version
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMPI"
# uncomment the next line if you want compressed tau-resolved data files
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DZLIB"
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMKL_DIRECT_CALL"

F90OPTFLAGS="-O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
# uncomment the next line if you want to use additional openmp parallelization
F90OPTFLAGS=${F90OPTFLAGS}" -parallel -qopenmp"
F90USEFULFLAGS="-cpp -std03"
export F90USEFULFLAGS
export F90OPTFLAGS

FL="${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
export FL

export Libs=${DIR}"/Libraries/"
#export LIB_BLAS_LAPACK="-lblas -llapack"
export LIB_BLAS_LAPACK="-mkl -lz"
#export LIB_BLAS_LAPACK="-mkl=parallel"
