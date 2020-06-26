STABCONFIGURATION=""
# STABCONFIGURATION=${STABCONFIGURATION}" -DQRREF"

PROGRAMMCONFIGURATION="-DMPI"
INTELCOMPILER="mpiifort"
GNUCOMPILER="mpifort"
MPICOMP=1

# default optimization flags for Intel compiler
INTELOPTFLAGS="-cpp -O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
INTELOPTFLAGS="-cpp -O3 "
# uncomment the next line if you want to use additional openmp parallelization
INTELOPTFLAGS="${INTELOPTFLAGS} -parallel -qopenmp"
INTELUSEFULFLAGS="-std03"

# default optimization flags for GNU compiler
GNUOPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
# uncomment the next line if you want to use additional openmp parallelization
GNUOPTFLAGS="${GNUOPTFLAGS} -fopenmp"
GNUUSEFULFLAGS="-std=f2003"

MACHINE=""
Machinev=0
MODE=""
modev=0
STAB=""
stabv=0

RED='\033[0;31m'
NC='\033[0m' # No Color

shopt -s nocasematch

while [ "$#" -gt "0" ]; do
  case "$1" in
    STAB1|STAB2|STAB3|LOG)
      if [ "$stabv" == "1" ]; then
         echo -e "Additional STAB configuration found. Overwriting "$STAB" with" $1 "."
      fi
      STAB="$1"
      stabv="1"
      shift 1
    ;;
    noMPI|MPI|Tempering|serial)
      if [ "$modev" == "1" ]; then
         echo -e "Additional MODE configuration found. Overwriting "$MODE" with" $1 "."
      fi
      MODE="$1"
      modev="1"
      shift 1
    ;;
    *)
      if [ "$Machinev" == "1" ]; then
         echo -e "Additional MACHINE / unrecognized configuration found. Overwriting "$MACHINE" with" $1 "."
      fi
      MACHINE="$1"
      Machinev="1"
      shift 1
    ;;  
  esac
done 


echo ""

case $MODE in

# STAB1|STAB2|STAB3|LOG)
#$3="$2"
#echo "Please provide MODE as second and STAB as third argument."
#echo "If you wish to keep the default mode, provide MPI as second argument"
#echo "Usage '. configureHPC.sh "$1" MODE "$2
#;;

noMPI|serial)
echo "serial job."
PROGRAMMCONFIGURATION=""
INTELCOMPILER="ifort"
GNUCOMPILER="gfortran"
MPICOMP=0
;;

Tempering)
echo "Activating parallel tempering."
echo "This requires also MPI parallization which is set as well."
PROGRAMMCONFIGURATION="-DMPI -DTEMPERING"
INTELCOMPILER="mpiifort"
GNUCOMPILER="mpifort"
MPICOMP=1
;;

MPI)
echo "Activating MPI parallization."
PROGRAMMCONFIGURATION="-DMPI"
INTELCOMPILER="mpiifort"
GNUCOMPILER="mpifort"
MPICOMP=1
;;

*)
echo -e "Activating "${RED}"MPI parallization (default)"${NC}"."
echo "To turn MPI off, pass noMPI as the second argument."
echo "To turn on parallel tempering, pass Tempering as the second argument."
PROGRAMMCONFIGURATION="-DMPI"
INTELCOMPILER="mpiifort"
GNUCOMPILER="mpifort"
MPICOMP=1
;;

esac

echo ""

case $STAB in

STAB1)
STABCONFIGURATION="${STABCONFIGURATION} -DSTAB1"
echo "Using older stabilization with UDV decompositions"
;;

STAB2)
STABCONFIGURATION="${STABCONFIGURATION} -DSTAB2"
echo "Using older stabilization with UDV decompositions and additional normalizations"
;;

STAB3)
STABCONFIGURATION="${STABCONFIGURATION} -DSTAB3"
echo "Using newest stabilization which seperates large and small scales"
;;

LOG)
STABCONFIGURATION="${STABCONFIGURATION} -DLOG"
echo "Using log storage for internal scales"
;;

# noMPI|MPI|Tempering)
#echo "Please provide MODE as second and STAB as third argument."
#echo "If you wish to keep the default mode, provide MPI as second argument"
#echo "Usage '. configureHPC.sh "$1" "$3" <STAB>"
#;;

*)
echo -e "Using "${RED}"default stabilization"${NC}
echo "Possible alternative options are STAB1, STAB2, STAB3 and LOG"
;;

esac

case $MACHINE in

#Fakhers MacBook
FakhersMAC)

# F90OPTFLAGS=$GNUOPTFLAGS
F90OPTFLAGS="$GNUOPTFLAGS -Wconversion -fcheck=all -g -fbacktrace"
F90USEFULFLAGS=$GNUUSEFULFLAGS 
if [ "$MPICOMP" -eq "0" ]; then
f90="gfortran"
else
f90="$mpif90"
fi
LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
;;

#Development
Devel|Development)

# F90OPTFLAGS="$GNUOPTFLAGS -Wconversion -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal"
F90OPTFLAGS=$GNUOPTFLAGS" -Wconversion -Werror -fcheck=all -g -fbacktrace "
# F90OPTFLAGS=$GNUOPTFLAGS" -Wconversion -Wcompare-reals -fcheck=all -g -fbacktrace "
F90USEFULFLAGS=$GNUUSEFULFLAGS

f90="gfortran"
LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
;;

#LRZ enviroment
SuperMUC)
module switch mpi.ibm  mpi.intel/2018
module switch intel intel/18.0
module switch mkl mkl/2018

F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
f90=mpiifort
LIB_BLAS_LAPACK=$MKL_LIB
;;

#LRZ enviroment
SuperMUC-NG|NG)
#module switch mpi.intel  mpi.intel/2018
#module switch intel intel/18.0
#module switch mkl mkl/2018
module load  mpi.intel
module load intel
module load mkl

F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
export f90=mpiifort
export LIB_BLAS_LAPACK=$MKL_LIB
;;

#JUWELS enviroment
JUWELS)
module load Intel
module load IntelMPI
module load imkl

F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
f90=mpiifort
LIB_BLAS_LAPACK="-mkl"
;;

#Intel (as Hybrid code)
Intel)
F90OPTFLAGS=$INTELOPTFLAGS
F90USEFULFLAGS=$INTELUSEFULFLAGS
f90=$INTELCOMPILER
LIB_BLAS_LAPACK="-mkl"
;;

#GNU (as Hybrid code)
GNU)
F90OPTFLAGS=$GNUOPTFLAGS
F90USEFULFLAGS=$GNUUSEFULFLAGS
f90=$GNUCOMPILER
LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
;;

#PGI
PGI)
if [ "$MPICOMP" -eq "0" ]; then
f90="pgfortran"
else
f90="$(dirname $(dirname $(command -v pgfortran)))/mpi/openmpi/bin/mpifort"
fi
LIB_BLAS_LAPACK="-llapack -lblas"
F90OPTFLAGS="-Mpreprocess -O3 -mp"
F90USEFULFLAGS="-Minform=inform"
;;

#Default (unknown machine)
*)
echo
echo -e ${RED}"   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"${NC}
echo -e ${RED}"   !!               UNKNOW MACHINE               !!"${NC}
echo -e ${RED}"   !!         IGNORING PARALLEL SETTINGS         !!"${NC}
echo -e ${RED}"   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"${NC}
echo 
echo "Activating fallback option with gfortran for SERIAL JOB."
echo 
echo "usage 'source configureHPC.sh MACHINE MODE STAB'"
echo 
echo "Please choose one of the following machines:"
echo " * SuperMUC"
echo " * SuperMUC-NG"
echo " * JUWELS"
echo " * Devel"
echo " * Intel"
echo " * GNU"
echo " * FakhersMAC"
echo "Possible modes are MPI (default), noMPI and Tempering"
echo "Possible stab are no-argument (default), STAB1 (old), STAB2 (old), STAB3 (newest)"
echo "    and LOG (increases accessible scales, e.g. in beta or interaction strength by solving NaN issues)"

PROGRAMMCONFIGURATION=""
F90OPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
F90USEFULFLAGS=""

f90=gfortran
LIB_BLAS_LAPACK="-llapack -lblas"
;;

esac

PROGRAMMCONFIGURATION="$STABCONFIGURATION $PROGRAMMCONFIGURATION"

Libs="$(pwd)/Libraries"
ALF_INC="-I${Libs}/Modules"
ALF_LIB="${Libs}/Modules/modules_90.a ${Libs}/libqrref/libqrref.a ${LIB_BLAS_LAPACK}"
export ALF_LIB

export ALF_DIR="$(pwd)"
export ALF_FC=$f90

if [ ! -z ${ALF_FLAGS_EXT+x} ]; then 
  echo; echo "Appending additional compiler flag '${ALF_FLAGS_EXT}'"
fi

ALF_FLAGS_QRREF="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
ALF_FLAGS_MODULES="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
ALF_FLAGS_ANA="${F90OPTFLAGS} ${ALF_INC} ${ALF_FLAGS_EXT}"
ALF_FLAGS_PROG="${F90USEFULFLAGS} ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_INC} ${ALF_FLAGS_EXT}"
export ALF_FLAGS_QRREF
export ALF_FLAGS_MODULES
export ALF_FLAGS_ANA
export ALF_FLAGS_PROG

echo
echo "To compile your program use:    'make TARGET'"
echo
