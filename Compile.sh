# -DMPI selects MPI.
# -DSTAB1    Singular value decomposition for stabilization
# -DSTAB2    QR  with pivoting. Packed form of QR factoriztion  is not used.
#  Default  stabilization QR with pivotting. Packed form of QR factoriztion  is used. 
PROGRAMMCONFIGURATION="-DMPI"
PROGRAMMCONFIGURATION=""
f90="gfortran"
export f90
F90OPTFLAGS="-O3 -Wconversion  -fcheck=all"
F90OPTFLAGS="-O3 "
export F90OPTFLAGS
F90USEFULFLAGS="-cpp -std=f2003"
F90USEFULFLAGS="-cpp "
export F90USEFULFLAGS
FL="-c ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
export FL
DIR=`pwd`
export DIR
Libs=${DIR}"/Libraries/"
export Libs
LIB_BLAS_LAPACK="-llapack -lblas"
export LIB_BLAS_LAPACK

cd  Libraries 
make 
cd ../Analysis
make
cd ../Prog
make Examples
