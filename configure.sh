# setting QRREF has the highest priority. Setting nothing selects System lapack for the QR decomposition.
# Setting OLDNAG selects syntax of NAG Versions before NAG Mark 17 (roughly....).
# In addition the NAG library has to be specified in LIB_BLAS_LAPACK.
# -DMPI selects MPI.
PROGRAMMCONFIGURATION="-DnoMPI"
DIR=`pwd`
export DIR
f90="gfortran"
export f90
FL="-c -O3 ${PROGRAMMCONFIGURATION}"
export FL
Libs=${DIR}"/Libraries/"
export Libs
LIB_BLAS_LAPACK="-llapack -lblas"
export LIB_BLAS_LAPACK
