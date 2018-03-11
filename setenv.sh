export DIR=`pwd`
f90=$mpif90
f90=gfortran
export f90
F90OPTFLAGS="-O3 -Wconversion  -fcheck=all"
F90OPTFLAGS="-O3 -Wconversion "
export F90OPTFLAGS
F90USEFULFLAGS=" -cpp -std=f2003"
F90USEFULFLAGS=" -cpp"
export F90USEFULFLAGS
FL="-c ${F90OPTFLAGS} ${PROGRAMCONFIGURATION}"
export FL
Libs="$DIR/Libraries/"
export Libs
LIB_BLAS_LAPACK=" -llapack -lblas"
export LIB_BLAS_LAPACK

