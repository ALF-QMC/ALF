export DIR=`pwd`
export f90="ifort"
export FL="-c  -w   -r8  -O3 "
export FL="-c  -w  -O3"
export FL="-c  -w   -fdefault-real-8 -O3   -fbounds-check"
export FL="-c  -w   -r8   -fbounds-check -ftrace=full -C -g"
export FL="-c  -w   -r8   -O3"
export FL="-c  -w   -r8  -O3"
export LIB_BLAS_LAPACK=${HOME}"/lib_90/LaPack/lapack.a  "${HOME}"/lib_90/Blas/libblas.a"
