# setting QRREF has the highest priority. Setting nothing selects System lapack for the QR decomposition.
# -DQRREF sets reference QR
# -DMPI selects MPI.
# -DSTAB1 selects an altrnative stabilitation scheme.
PROGRAMMCONFIGURATION=""
f90="gfortran"
export f90
F90OPTFLAGS="-O3"
export F90OPTFLAGS
FL="-c ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
export FL
DIR=`pwd`
export DIR
Libs=${DIR}"/Libraries/"
export Libs
LIB_BLAS_LAPACK="-llapack -lblas"
export LIB_BLAS_LAPACK
