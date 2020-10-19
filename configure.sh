#!/bin/sh
STABCONFIGURATION=""
# STABCONFIGURATION="${STABCONFIGURATION} -DQRREF"

# default optimization flags for Intel compiler
INTELOPTFLAGS="-cpp -O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
INTELOPTFLAGS="-cpp -O3"
INTELOPTFLAGS="$INTELOPTFLAGS -no-wrap-margin"
#INTELOPTFLAGS="$INTELOPTFLAGS -traceback"
# uncomment the next line if you want to use additional openmp parallelization
INTELOPTFLAGS="${INTELOPTFLAGS} -parallel -qopenmp"
INTELUSEFULFLAGS="-std08"

# default optimization flags for GNU compiler
GNUOPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math -fmax-errors=10"
# uncomment the next line if you want to use additional openmp parallelization
GNUOPTFLAGS="${GNUOPTFLAGS} -fopenmp"
GNUUSEFULFLAGS="-std=f2008"

MACHINE=""
Machinev=0
MODE=""
modev=0
STAB=""
stabv=0

RED='\033[0;31m'
NC='\033[0m' # No Color

while [ "$#" -gt "0" ]; do
  ARG="$(echo "$1" | tr '[:lower:]' '[:upper:]')"
  shift 1
  case "$ARG" in
    STAB1|STAB2|STAB3|LOG)
      if [ "$stabv" = "1" ]; then
         printf "Additional STAB configuration found. Overwriting %s with %s .\n" "$STAB" "$ARG"
      fi
      STAB="$ARG"
      stabv="1"
    ;;
    NOMPI|MPI|TEMPERING|SERIAL)
      if [ "$modev" = "1" ]; then
         printf "Additional MODE configuration found. Overwriting %s with %s .\n" "$MODE" "$ARG"
      fi
      MODE="$ARG"
      modev="1"
    ;;
    *)
      if [ "$Machinev" = "1" ]; then
         printf "Additional MACHINE / unrecognized configuration found. Overwriting %s with %s .\n" "$MACHINE" "$ARG"
      fi
      MACHINE="$ARG"
      Machinev="1"
    ;;
  esac
done

printf "\n"

case $MODE in
  NOMPI|SERIAL)
    printf "serial job.\n"
    PROGRAMMCONFIGURATION=""
    INTELCOMPILER="ifort"
    GNUCOMPILER="gfortran"
    MPICOMP=0
  ;;

  TEMPERING)
    printf "Activating parallel tempering.\n"
    printf "This requires also MPI parallization which is set as well.\n"
    PROGRAMMCONFIGURATION="-DMPI -DTEMPERING"
    INTELCOMPILER="mpiifort"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;

  MPI)
    printf "Activating MPI parallization.\n"
    PROGRAMMCONFIGURATION="-DMPI"
    INTELCOMPILER="mpiifort"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;

  *)
    printf "Activating ${RED}MPI parallization (default)${NC}.\n"
    printf "To turn MPI off, pass noMPI as the second argument.\n"
    printf "To turn on parallel tempering, pass Tempering as the second argument.\n"
    PROGRAMMCONFIGURATION="-DMPI"
    INTELCOMPILER="mpiifort"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;
esac

printf "\n"

case $STAB in
  STAB1)
    STABCONFIGURATION="${STABCONFIGURATION} -DSTAB1"
    printf "Using older stabilization with UDV decompositions\n"
  ;;

  STAB2)
    STABCONFIGURATION="${STABCONFIGURATION} -DSTAB2"
    printf "Using older stabilization with UDV decompositions and additional normalizations\n"
  ;;

  STAB3)
    STABCONFIGURATION="${STABCONFIGURATION} -DSTAB3"
    printf "Using newest stabilization which seperates large and small scales\n"
  ;;

  LOG)
    STABCONFIGURATION="${STABCONFIGURATION} -DLOG"
    printf "Using log storage for internal scales\n"
  ;;

  *)
    printf "Using ${RED}default stabilization${NC}\n"
    printf "Possible alternative options are STAB1, STAB2, STAB3 and LOG\n"
  ;;
esac

case $MACHINE in
  #Fakhers MacBook
  FAKHERSMAC)
    # F90OPTFLAGS=$GNUOPTFLAGS
    F90OPTFLAGS="$GNUOPTFLAGS -Wconversion  -Wuninitialized  -fcheck=all -g -fbacktrace"
    F90USEFULFLAGS="$GNUUSEFULFLAGS"
    if [ "$MPICOMP" -eq "0" ]; then
    ALF_FC="gfortran"
    else
    ALF_FC="$mpif90"
    fi
    LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
  ;;

  #Development
  DEVEL|DEVELOPMENT)
    # F90OPTFLAGS="$GNUOPTFLAGS -Wconversion -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal"
    F90OPTFLAGS="$GNUOPTFLAGS -Wconversion -Werror -fcheck=all -g -fbacktrace "
    # F90OPTFLAGS=$GNUOPTFLAGS" -Wconversion -Wcompare-reals -fcheck=all -g -fbacktrace "
    F90USEFULFLAGS="$GNUUSEFULFLAGS"

    ALF_FC="$GNUCOMPILER"
    LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
  ;;

  #LRZ enviroment
  SUPERMUC)
    module switch mpi.ibm  mpi.intel/2018
    module switch intel intel/18.0
    module switch mkl mkl/2018

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="$MKL_LIB"
  ;;

  #LRZ enviroment
  SUPERMUC-NG|NG)
    #module switch mpi.intel  mpi.intel/2018
    #module switch intel intel/18.0
    #module switch mkl mkl/2018
    module load  mpi.intel
    module load intel
    module load mkl

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="$MKL_LIB"
  ;;

  #JUWELS enviroment
  JUWELS)
    module load Intel
    module load IntelMPI
    module load imkl

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="-mkl"
  ;;

  #Intel (as Hybrid code)
  INTEL)
    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="$INTELCOMPILER"
    LIB_BLAS_LAPACK="-mkl"
  ;;

  #GNU (as Hybrid code)
  GNU)
    F90OPTFLAGS="$GNUOPTFLAGS"
    F90USEFULFLAGS="$GNUUSEFULFLAGS"
    ALF_FC="$GNUCOMPILER"
    LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
  ;;

  #PGI
  PGI)
    if [ "$MPICOMP" -eq "0" ]; then
      ALF_FC="pgfortran"
    else
      ALF_FC="mpifort"
      printf "\n${RED}   !! Compiler set to 'mpifort' !!\n"
      printf "If this is not your PGI MPI compiler you have to set it manually through:\n"
      printf "    'export ALF_FC=<mpicompiler>'${NC}\n"
    fi
    LIB_BLAS_LAPACK="-llapack -lblas"
    F90OPTFLAGS="-Mpreprocess -O1 -mp"
    F90USEFULFLAGS="-Minform=inform"
  ;;

  #Default (unknown machine)
  *)
    printf "\n"
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n"
    printf "${RED}   !!               UNKNOW MACHINE               !!${NC}\n"
    printf "${RED}   !!         IGNORING PARALLEL SETTINGS         !!${NC}\n"
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n"
    printf "\n"
    printf "Activating fallback option with gfortran for SERIAL JOB.\n"
    printf "\n"
    printf "usage 'source configureHPC.sh MACHINE MODE STAB'\n"
    printf "\n"
    printf "Please choose one of the following machines:\n"
    printf " * SuperMUC\n"
    printf " * SuperMUC-NG\n"
    printf " * JUWELS\n"
    printf " * Devel\n"
    printf " * Intel\n"
    printf " * GNU\n"
    printf " * FakhersMAC\n"
    printf "Possible modes are MPI (default), noMPI and Tempering\n"
    printf "Possible stab are no-argument (default), STAB1 (old), STAB2 (old), STAB3 (newest)\n"
    printf "and LOG (increases accessible scales, e.g. in beta or interaction strength by solving NaN issues)\n"

    PROGRAMMCONFIGURATION=""
    F90OPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
    F90USEFULFLAGS=""

    ALF_FC="gfortran"
    LIB_BLAS_LAPACK="-llapack -lblas"
  ;;
esac

PROGRAMMCONFIGURATION="$STABCONFIGURATION $PROGRAMMCONFIGURATION"

Libs="$PWD/Libraries"
ALF_INC="-I${Libs}/Modules"
ALF_LIB="${Libs}/Modules/modules_90.a ${Libs}/libqrref/libqrref.a ${LIB_BLAS_LAPACK}"
export ALF_LIB

export ALF_DIR="$PWD"
export ALF_FC

if [ ! -z "${ALF_FLAGS_EXT+x}" ]; then
  printf "\nAppending additional compiler flag '%s'\n" "${ALF_FLAGS_EXT}"
fi

ALF_FLAGS_QRREF="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
ALF_FLAGS_MODULES="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
ALF_FLAGS_ANA="${F90USEFULFLAGS} ${F90OPTFLAGS} ${ALF_INC} ${ALF_FLAGS_EXT}"
ALF_FLAGS_PROG="${F90USEFULFLAGS} ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_INC} ${ALF_FLAGS_EXT}"
export ALF_FLAGS_QRREF
export ALF_FLAGS_MODULES
export ALF_FLAGS_ANA
export ALF_FLAGS_PROG

printf "\nTo compile your program use:    'make TARGET'\n\n"
