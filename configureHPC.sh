PROGRAMMCONFIGURATION=""
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DSTAB1"
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DSTAB2"
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DQRREF"

# default optimization flags for Intel compiler
F90OPTFLAGS="-O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
# uncomment the next line if you want to use additional openmp parallelization
F90OPTFLAGS=${F90OPTFLAGS}" -parallel -qopenmp"
F90USEFULFLAGS="-cpp"

export DIR=`pwd`

case $2 in

noMPI)
echo "seriell job."
;;

Tempering)
echo "Activating parallel tempering."
echo "This requires also MPI parallization which is set as well."
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMPI -DTEMPERING"
;;

MPI|*)
echo "Activating MPI parallization (default)."
echo "To turn MPI off, pass noMPI as the second argument."
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMPI"
;;

esac

echo ""

case $1 in

#Development
Devel)

PROGRAMMCONFIGURATION=""
F90OPTFLAGS="-O3 -ffree-line-length-none -Wconversion -Werror"
F90USEFULFLAGS="-cpp"

export f90=gfortran
export LIB_BLAS_LAPACK="-llapack -lblas"
;;

#LRZ enviroment
SuperMUC)
module switch mpi.ibm mpi.intel
module switch intel intel/17.0
module switch mkl mkl/2017

export f90=mpif90
export LIB_BLAS_LAPACK=$MKL_LIB
;;

#JURECA enviroment
JURECA)
module load Intel
module load IntelMPI
module load imkl

export f90=mpiifort
export LIB_BLAS_LAPACK="-mkl"
;;

#Intel (as Hybrid code)
Intel)
export f90=mpiifort
export LIB_BLAS_LAPACK="-mkl"
;;

#Matrix23 PGI
Matrix23)
export f90=pgfortran
export LIB_BLAS_LAPACK="-L/opt/pgi/linux86-64/17.4/lib -llapack -lblas"
F90OPTFLAGS="-O3 -mp"
F90USEFULFLAGS="-Mpreprocess -Minform=inform"
;;

#Default (unknown machine)
*)
echo "Please choose one of the following machines:"
echo "SuperMUC"
echo "JURECA"
echo
echo "usage 'source configureHPC.sh MACHINE'"
echo 
echo "Activating fallback option with gfortran for SERIAL JOB."

PROGRAMMCONFIGURATION=""
F90OPTFLAGS="-O3 -ffree-line-length-none"
F90OPTFLAGS="-O3 -ffree-line-length-none  -fcheck=all"
F90OPTFLAGS="-O3 -ffree-line-length-none  -Wconversion -fcheck=all"
F90OPTFLAGS="-O3 -ffree-line-length-none  -Wconversion "
F90USEFULFLAGS="-cpp"

export f90="gfortran"
export LIB_BLAS_LAPACK="-llapack -lblas"
;;

esac

export F90USEFULFLAGS
export F90OPTFLAGS

FL="-c ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
export FL

export Libs=${DIR}"/Libraries/"

echo
echo "To compile your program use:    'make -f MakefileHPC TARGET'"



