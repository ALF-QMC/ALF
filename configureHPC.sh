PROGRAMMCONFIGURATION=""
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DSTAB1"
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DSTAB2"
# PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DQRREF"

# uncomment the next line if you want an parallel tempering version
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DTEMPERING"

# uncomment the next line if you want an MPI parallel version
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMPI"

F90OPTFLAGS="-O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
# uncomment the next line if you want to use additional openmp parallelization
F90OPTFLAGS=${F90OPTFLAGS}" -parallel -qopenmp"
F90USEFULFLAGS="-cpp"

export DIR=`pwd`

case $1 in

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

#Default (unknown machine)
*)
echo "Please choose one of the following machines:"
echo "SuperMUC"
echo "JURECA"
echo
echo "usage 'source configureHPC.sh MACHINE'"
echo 
echo "Activating fallback option with gfortran for SERIAL JOB."

#PROGRAMMCONFIGURATION=""
F90OPTFLAGS="-O3 -ffree-line-length-none  -fcheck=all"
F90OPTFLAGS="-O3 -ffree-line-length-none  -Wconversion"
F90USEFULFLAGS="-cpp"

export f90=$mpif90
#"gfortran"
# "gfortran"
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



