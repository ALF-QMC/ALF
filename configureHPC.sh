PROGRAMMCONFIGURATION=""
STABCONFIGURATION=""
# STABCONFIGURATION=${STABCONFIGURATION}" -DSTAB1"
# STABCONFIGURATION=${STABCONFIGURATION}" -DSTAB2"
# STABCONFIGURATION=${STABCONFIGURATION}" -DQRREF"

# default optimization flags for Intel compiler
INTELOPTFLAGS="-O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
# uncomment the next line if you want to use additional openmp parallelization
INTELOPTFLAGS=${INTELOPTFLAGS}" -parallel -qopenmp"
INTELUSEFULFLAGS="-cpp -std03"

# default optimization flags for GNU compiler
GNUOPTFLAGS="-O3 -ffree-line-length-none -ffast-math"
# uncomment the next line if you want to use additional openmp parallelization
GNUOPTFLAGS=${GNUOPTFLAGS}" -fopenmp"
GNUUSEFULFLAGS="-cpp -std=f2003"

export DIR=`pwd`

case $2 in

noMPI)
echo "seriell job."
INTELCOMPILER="ifort"
GNUCOMPILER="gfortran"
FAKHERCOMPILER="gfortran"
;;

Tempering)
echo "Activating parallel tempering."
echo "This requires also MPI parallization which is set as well."
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMPI -DTEMPERING"
INTELCOMPILER="mpiifort"
# INTELUSEFULFLAGS="-cpp"
GNUCOMPILER="mpifort"
# GNUUSEFULFLAGS="-cpp"
FAKHERCOMPILER=$(mpif90)
;;

MPI|*)
echo "Activating MPI parallization (default)."
echo "To turn MPI off, pass noMPI as the second argument."
echo "To turn use parallel tempering, pass Tempering as the second argument."
PROGRAMMCONFIGURATION=${PROGRAMMCONFIGURATION}" -DMPI"
INTELCOMPILER="mpiifort"
# INTELUSEFULFLAGS="-cpp"
GNUCOMPILER="mpifort"
# GNUUSEFULFLAGS="-cpp"
FAKHERCOMPILER=$(mpif90)
;;

esac

echo ""

case $1 in

#Fakhers MacBook
FakhersMAC)

F90OPTFLAGS=$GNUOPTFLAGS
F90USEFULFLAGS=$GNUUSEFULFLAGS

export f90=$FAKHERCOMPILER
export LIB_BLAS_LAPACK="-llapack -lblas"
;;

#Development
Devel)

F90OPTFLAGS=$GNUOPTFLAGS" -Wconversion -Werror -fcheck=all"
F90USEFULFLAGS=$GNUUSEFULFLAGS

export f90=$GNUCOMPILER
export LIB_BLAS_LAPACK="-llapack -lblas"
;;

#LRZ enviroment
SuperMUC)
module switch mpi.ibm mpi.intel
module switch intel intel/17.0
module switch mkl mkl/2017

F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
export f90=mpif90
export LIB_BLAS_LAPACK=$MKL_LIB
;;

#JURECA enviroment
JURECA)
module load Intel
module load IntelMPI
module load imkl

F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
export f90=mpiifort
export LIB_BLAS_LAPACK="-mkl"
;;

#Intel (as Hybrid code)
Intel)
F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
export f90=$INTELCOMPILER
export LIB_BLAS_LAPACK="-mkl"
;;

#GNU (as Hybrid code)
GNU)
F90OPTFLAGS=$GNUOPTFLAGS
F90USEFULFLAGS=$GNUUSEFULFLAGS
export f90=$GNUCOMPILER
export LIB_BLAS_LAPACK="-llapack -lblas"
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
echo "usage 'source configureHPC.sh MACHINE MODE'"
echo 
echo "Please choose one of the following machines:"
echo "  SuperMUC"
echo "  JURECA"
echo "  Devel"
echo "  Intel"
echo "  GNU"
echo
echo "Possible modes are MPI (default), noMPI and Tempering"
echo
echo "Please choose one of the following machines:"
echo "Activating fallback option with gfortran for SERIAL JOB."

PROGRAMMCONFIGURATION=""
F90OPTFLAGS="-O3 -ffree-line-length-none -ffast-math"
F90USEFULFLAGS="-cpp"

export f90=gfortran
export LIB_BLAS_LAPACK="-llapack -lblas"
;;

esac

PROGRAMMCONFIGURATION=$STABCONFIGURATION" "$PROGRAMMCONFIGURATION

export F90USEFULFLAGS
export F90OPTFLAGS

FL="-c ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
export FL

export Libs=${DIR}"/Libraries/"

echo
echo "To compile your program use:    'make TARGET'"



