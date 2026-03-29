#!/bin/sh
# shellcheck disable=SC2059
# This script sets necessary environment variables for compiling ALF.
# You need to source it prior to executing make.
USAGE="usage 'source configure.sh MACHINE MODE STAB'

Please choose one of the following MACHINEs:
 * GNU
 * Intel
 * IntelLLVM or IntelX
 * PGI
 * SuperMUC-NG
 * JUWELS
 * FRITZ
Possible MODEs are:
 * MPI (default)
 * noMPI
 * Tempering
 * PARALLEL_PARAMS (shorthand PP)
Possible STABs are:
 * <no-argument> (default)
 * STAB1 (old)
 * STAB2 (old)
 * STAB3 (newest)
 * LOG (increases accessible scales, e.g. in beta or interaction strength by solving NaN issues)
Further optional arguments: 
  Devel: Compile with additional flags for development and debugging
  HDF5: Compile with HDF5
  NO-INTERACTIVE: Do not ask for user confirmation during excution of this script
  NO-FALLBACK: Do not use a fallback option in case of an unknown/no machine,
               but instead return with value 1
To hand an additional flag to the compiler, export it in the varible ALF_FLAGS_EXT prior to sourcing this script.

ALF usually self-compiles HDF5 and stores the library in subdirectories of ALF/HDF5.
This behavior can be changed by setting the environment variable ALF_HDF5_DIR.

For more details check the documentation.\n"

STABCONFIGURATION=""
# STABCONFIGURATION="${STABCONFIGURATION} -DQRREF"

export ALF_DIR="$PWD"

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Create temporary directory for various checks with temporary files to be run in parallel
tmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir')
printf "\n${GREEN}Temporary directory %s created${NC}\n" "$tmpdir"

set_hdf5_flags()
{
  CC="$1" FC="$2" CXX="$3"
  
  $FC -o "$tmpdir/get_compiler_version.out" get_compiler_version.F90
  compiler_vers=$("$tmpdir/get_compiler_version.out" | sed 's/[ ,()]/_/g')
  
  H5_major=1
  H5_minor=14
  H5_patch=5
  H5_suff=""
  if [ -n "${ALF_HDF5_DIR+x}" ]; then
    printf "\nUsing custom HDF5 directory '%s'\n" "${ALF_HDF5_DIR}"
    HDF5_DIR="${ALF_HDF5_DIR}/${compiler_vers}"
  else
    HDF5_DIR="$ALF_DIR/HDF5/${compiler_vers}"
  fi
  if [ ! -d "$HDF5_DIR" ]; then
    printf "\nHDF5 is not yet installed for compiler '%s'.\n" "$compiler_vers"
    printf "ALF does never use global HDF5 libraries, but installs it locally in subfolders of '%s/HDF5'.\n" "$ALF_DIR"
    if [ "$NO_INTERACTIVE" = "" ]; then
      printf "Do you want download and install it now locally in the ALF folder? (Y/n):"
      read -r yn
    else
      yn="Y"
    fi
    case "$yn" in
      y|Y|"")
        printf "${GREEN}Downloading and installing HDF5 in %s.${NC}\n" "$HDF5_DIR"
        CC="$CC" FC="$FC" CXX="$CXX" HDF5_DIR="$HDF5_DIR" "$ALF_DIR/HDF5/install_hdf5.sh" ${H5_major} ${H5_minor} ${H5_patch} "${H5_suff}" || return 1
      ;;
      *) 
        printf "Skipping installation of HDF5.\n"
        rm -r "$tmpdir"
        printf "\n${GREEN}Temporary directory %s deleted${NC}\n" "$tmpdir"
        return 1
      ;;
    esac
  fi
  INC_HDF5="-I$HDF5_DIR/include"
  LIB_HDF5="-L$HDF5_DIR/lib $HDF5_DIR/lib/libhdf5hl_fortran.a $HDF5_DIR/lib/libhdf5_hl.a"
  EXTRA_LIBRARIES="$("$HDF5_DIR"/bin/h5fc -showconfig | grep 'Extra libraries:' | cut -f2 -d':')"
  # EXTRA_LIBRARIES="-lz -ldl -lm"
  LIB_HDF5="$LIB_HDF5 $HDF5_DIR/lib/libhdf5_fortran.a $HDF5_DIR/lib/libhdf5.a $EXTRA_LIBRARIES -Wl,-rpath -Wl,$HDF5_DIR/lib"

  if ! "$HDF5_DIR/bin/h5fc" -showconfig | grep "deflate(zlib)" > /dev/null; then
    printf "${RED}Warning: HDF5 installed without compression capabilies. The output files will not be compressed!${NC}\n" 1>&2
  fi
}

check_libs()
{
    FC="$1" LIBS="$2"
    FC0="$(echo "$FC" | cut -f1 -d' ')"
    if command -v "$FC0" > /dev/null; then       # Compiler binary found
        if sh -c "$FC check_libs.f90 $LIBS -o $tmpdir/check_libs.out"; then  # Compiling with $LIBS is successful
            "$tmpdir/check_libs.out" || (
              printf "${RED}\n==== Error: Execution of test program using compiler <%s> ====${NC}\n" "$FC" 1>&2
              printf "${RED}==== and linear algebra libraries <%s> not successful. ====${NC}\n\n" "$LIBS" 1>&2
              # script gets terminated, so remove tmpdir
              rm -r "$tmpdir"
              printf "\n${GREEN}Temporary directory %s deleted${NC}\n" "$tmpdir"
              return 1
              )
        else
            printf "${RED}\n==== Error: Linear algebra libraries <%s> not found. ====${NC}\n\n" "$LIBS" 1>&2
              # script gets terminated, so remove tmpdir
              rm -r "$tmpdir"
              printf "\n${GREEN}Temporary directory %s deleted${NC}\n" "$tmpdir"
            return 1
        fi
    else
        printf "${RED}\n==== Error: Compiler <%s> not found. ====${NC}\n\n" "$FC" 1>&2
        # script gets terminated, so remove tmpdir
        rm -r "$tmpdir"
        printf "\n${GREEN}Temporary directory %s deleted${NC}\n" "$tmpdir"
        return 1
    fi
}

check_python()
{
    if ! command -v python3 > /dev/null; then
        printf "${RED}\n==== Error: Python 3 not found. =====${NC}\n\n" 1>&2
        return 1
    fi
}

find_mkl_flag()
{
  if command -v ifort > /dev/null; then
    # default optimization flags for Intel compiler
    ifort_major=2000
    ifort_minor=0
    ifort_major=$(ifort --version | head -n 1 | awk '{print $(NF - 1)}' | cut -d '.' -f 1)
    ifort_minor=$(ifort --version | head -n 1 | awk '{print $(NF - 1)}' | cut -d '.' -f 2)
    if [ "$ifort_major" -gt 2021 ] || { [ "$ifort_major" -eq 2021 ] && [ "$ifort_minor" -gt 3 ]; }; then
      INTELMKL="-qmkl"
    else
      INTELMKL="-mkl"
    fi
  elif command -v ifx > /dev/null; then
    INTELMKL="-qmkl"
  else 
    printf "${RED}\n==== Error: MKL only supported for ifort compiler. ====${NC}\n\n" "$FC" 1>&2
  fi
}

set_intelcc()
{
  if command -v icx > /dev/null; then
    INTELCC="icx"
  elif command -v icc > /dev/null; then
    INTELCC="icc"
  elif command -v gcc > /dev/null; then
    INTELCC="gcc"
  else
    printf "${RED}\n==== Error: C compiler needed for HDF5. None of 'icx', 'icc', 'gcc' found ====${NC}\n\n" 1>&2
  fi
}

set_intelcxx()
{
  if command -v icpx > /dev/null; then
    INTELCXX="icpx"
  elif command -v icpc > /dev/null; then
    INTELCXX="icpc"
  elif command -v g++ > /dev/null; then
    INTELCXX="g++"
  else
    printf "${RED}\n==== Error: C++ compiler needed for HDF5. None of 'icpx', 'icpc', 'g++' found ====${NC}\n\n" 1>&2
  fi
}

# default optimization flags for Intel compiler
INTELOPTFLAGS="-cpp -O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
# INTELOPTFLAGS="-cpp -O3 "
#INTELOPTFLAGS="$INTELOPTFLAGS -traceback"
# uncomment the next line if you want to use additional openmp parallelization
INTELOPTFLAGS="${INTELOPTFLAGS} -parallel -qopenmp"
INTELDEVFLAGS="-warn all -check all -g -traceback"
INTELUSEFULFLAGS="-std08"

INTELLLVMOPTFLAGS="-cpp -O3"
INTELLLVMOPTFLAGS="-cpp -O3 -fp-model=fast=2 -no-prec-div -static -xHost -unroll -finline-functions -no-wrap-margin"
# uncomment the next line if you want to use additional openmp parallelization
# INTELLLVMOPTFLAGS="${INTELLLVMOPTFLAGS} -qopenmp"
INTELLLVMDEVFLAGS="-warn all -check all,nouninit -g -traceback"
INTELLLVMUSEFULFLAGS="-std08"


# default optimization flags for GNU compiler
GNUOPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
#GNUOPTFLAGS="-cpp -O0 -ffree-line-length-none"
# uncomment the next line if you want to use additional openmp parallelization
GNUOPTFLAGS="${GNUOPTFLAGS} -fopenmp"
# GNUDEVFLAGS="-Wconversion -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal"
GNUDEVFLAGS="-Wconversion -Werror -Wno-error=cpp -fcheck=all -g -fbacktrace -fmax-errors=10"
GNUUSEFULFLAGS="-std=f2008"

# default optimization flags for PGI compiler
PGIOPTFLAGS="-Mpreprocess -O3 -Mfprelaxed -fast"
# uncomment the next line if you want to use additional openmp parallelization
PGIOPTFLAGS="${PGIOPTFLAGS} -mp"
PGIDEVFLAGS="-Minform=inform -C -g -traceback"
PGIUSEFULFLAGS=""

MACHINE=""
Machinev=0
MODE=""
modev=0
STAB=""
stabv=0
HDF5_ENABLED=""
NO_INTERACTIVE=""
NO_FALLBACK=""

while [ "$#" -gt "0" ]; do
  ARG="$(echo "$1" | tr '[:lower:]' '[:upper:]')"
  shift 1
  case "$ARG" in
    STAB1|STAB2|STAB3|LOG)
      if [ "$stabv" = "1" ]; then
         printf "Additional STAB configuration found. Overwriting %s with %s .\n" "$STAB" "$ARG" 1>&2
      fi
      STAB="$ARG"
      stabv="1"
    ;;
    NOMPI|MPI|TEMPERING|SERIAL|PARALLEL_PARAMS|PP)
      if [ "$modev" = "1" ]; then
         printf "Additional MODE configuration found. Overwriting %s with %s .\n" "$MODE" "$ARG" 1>&2
      fi
      MODE="$ARG"
      modev="1"
    ;;
    HDF5)
      HDF5_ENABLED="1"
    ;;
    DEVEL|DEVELOPMENT)
      #DEVEL="1"
      GNUOPTFLAGS="$GNUOPTFLAGS $GNUDEVFLAGS"
      INTELOPTFLAGS="$INTELOPTFLAGS $INTELDEVFLAGS"
      INTELLLVMOPTFLAGS="$INTELLLVMOPTFLAGS $INTELLLVMDEVFLAGS"
      PGIOPTFLAGS="$PGIOPTFLAGS $PGIDEVFLAGS"
    ;;
    NO-INTERACTIVE)
      NO_INTERACTIVE="1"
    ;;
    NO-FALLBACK)
      NO_FALLBACK="1"
    ;;
    *)
      if [ "$Machinev" = "1" ]; then
         printf "Additional MACHINE / unrecognized configuration found. Overwriting %s with %s .\n" "$MACHINE" "$ARG" 1>&2
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
    INTELLLVMCOMPILER="ifx"
    GNUCOMPILER="gfortran"
    MPICOMP=0
  ;;

  TEMPERING)
    printf "Activating parallel tempering.\n"
    printf "This requires also MPI parallization which is set as well.\n"
    PROGRAMMCONFIGURATION="-DMPI -DTEMPERING"
    INTELCOMPILER="mpiifort"
    INTELLLVMCOMPILER="mpiifort -fc=ifx"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;

  MPI)
    printf "Activating MPI parallization.\n"
    PROGRAMMCONFIGURATION="-DMPI"
    INTELCOMPILER="mpiifort"
    INTELLLVMCOMPILER="mpiifort -fc=ifx"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;
 
  PARALLEL_PARAMS|PP)
    printf "Activating parallel runs with different parameters.\n"
    printf "This requires also MPI parallization which is set as well.\n"
    PROGRAMMCONFIGURATION="-DMPI -DTEMPERING -DPARALLEL_PARAMS"
    INTELCOMPILER="mpiifort"
    INTELLLVMCOMPILER="mpiifort -fc=ifx"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;

  *)
    printf "Activating ${RED}MPI parallization (default)${NC}.\n"
    printf "To turn MPI off, pass noMPI as the second argument.\n"
    printf "To turn on parallel tempering, pass Tempering as the second argument.\n"
    PROGRAMMCONFIGURATION="-DMPI"
    INTELCOMPILER="mpiifort"
    INTELLLVMCOMPILER="mpiifort -fc=ifx"
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
    STABCONFIGURATION="${STABCONFIGURATION} -DSTABLOG"
    printf "Using log storage for internal scales\n"
  ;;

  *)
    printf "Using ${RED}default stabilization${NC}\n"
    printf "Possible alternative options are STAB1, STAB2, STAB3 and LOG\n"
  ;;
esac

case $MACHINE in
  #GNU (as Hybrid code)
  GNU)
    F90OPTFLAGS="$GNUOPTFLAGS"
    F90USEFULFLAGS="$GNUUSEFULFLAGS"
    ALF_FC="$GNUCOMPILER"
    LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags gcc gfortran g++ || return 1
    fi
  ;;

  #Intel (as Hybrid code)
  INTEL)
    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="$INTELCOMPILER"
    find_mkl_flag || return 1
    LIB_BLAS_LAPACK="${INTELMKL}"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_intelcc
      set_intelcxx
      set_hdf5_flags "$INTELCC" ifort "$INTELCXX" || return 1
    fi
  ;;

  #Intel (as Hybrid code)
  INTELLLVM|INTELX)
    F90OPTFLAGS="$INTELLLVMOPTFLAGS"
    F90USEFULFLAGS="$INTELLLVMUSEFULFLAGS"
    ALF_FC="$INTELLLVMCOMPILER"
    INTELMKL="-qmkl"
    LIB_BLAS_LAPACK="${INTELMKL}"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_intelcc
      set_intelcxx
      set_hdf5_flags "$INTELCC" ifx "$INTELCXX" || return 1
    fi
  ;;

  #PGI
  PGI)
    F90OPTFLAGS="$PGIOPTFLAGS"
    F90USEFULFLAGS="$PGIUSEFULFLAGS"
    if [ "$MPICOMP" -eq "0" ]; then
      ALF_FC="pgfortran"
    else
      ALF_FC="mpifort"
      printf "\n${RED}   !! Compiler set to 'mpifort' !!\n" 1>&2
      printf "If this is not your PGI MPI compiler you have to set it manually through e.g.\n" 1>&2
      printf "    'export ALF_FC=<mpicompiler>'${NC}\n" 1>&2
    fi
    LIB_BLAS_LAPACK="-llapack -lblas"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags pgcc pgfortran pgc++ || return 1
    fi

  ;;

  #LRZ enviroment
  SUPERMUC-NG|NG)
    module load hdf5/1.10.7-intel21
    printf "\n${RED}   !!   unsetting  FORT_BLOCKSIZE  !!${NC}\n" 1>&2
    unset FORT_BLOCKSIZE

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="$MKL_LIB"
    LIB_HDF5="$HDF5_F90_SHLIB $HDF5_SHLIB"
    INC_HDF5="$HDF5_INC"
  ;;

  #JUWELS enviroment
  JUWELS)
    module load Intel
    module load IntelMPI
    module load imkl
    module load HDF5/1.10.6

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    find_mkl_flag || return 1
    LIB_BLAS_LAPACK="${INTELMKL}"
    LIB_HDF5="–lh5df_fortran"
    INC_HDF5=""
  ;;

  #NHR@FAU Fritz cluster
  FRITZ)
    module load intel
    module load intelmpi
    module load mkl

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="$INTELCOMPILER"
    find_mkl_flag || return 1
    LIB_BLAS_LAPACK="${INTELMKL}"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_intelcc
      set_intelcxx
      set_hdf5_flags "$INTELCC" ifort "$INTELCXX" || return 1
    fi
  ;;
  #Default (unknown machine)
  *)
    if [ "$NO_FALLBACK" = "1" ]; then
      printf "${RED}  !!     UNKNOW MACHINE     !!${NC}\n" 1>&2
      printf "${RED}  !!  exiting configure.sh  !!${NC}\n" 1>&2
      return 1
    fi
    printf "\n" 1>&2
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n" 1>&2
    printf "${RED}   !!               UNKNOW MACHINE               !!${NC}\n" 1>&2
    printf "${RED}   !!         IGNORING PARALLEL SETTINGS         !!${NC}\n" 1>&2
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n" 1>&2
    printf "\n" 1>&2
    printf "Activating fallback option with gfortran for SERIAL JOB - Deactivating MPI.\n" 1>&2
    printf "\n" 1>&2
    printf "$USAGE"
    PROGRAMMCONFIGURATION=""
    F90OPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
    F90USEFULFLAGS=""

    ALF_FC="gfortran"
    LIB_BLAS_LAPACK="-llapack -lblas"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags gcc gfortran g++ || return 1
    fi
  ;;
esac

check_libs "$ALF_FC" "${LIB_BLAS_LAPACK}" || return 1

check_python || return 1

PROGRAMMCONFIGURATION="$STABCONFIGURATION $PROGRAMMCONFIGURATION"

Libs="$ALF_DIR/Libraries"
ALF_INC="-I${Libs}/Modules"
ALF_LIB="${Libs}/Modules/modules_90.a ${LIB_BLAS_LAPACK} ${Libs}/libqrref/libqrref.a"
if [ "${HDF5_ENABLED}" = "1" ]; then
  echo; echo "HDF5 enabled"
  ALF_INC="${ALF_INC} ${INC_HDF5}"
  ALF_LIB="${ALF_LIB} ${LIB_HDF5}"
else
  echo; echo "HDF5 disabled"
fi
export ALF_LIB

export ALF_DIR
export ALF_FC="$ALF_FC"

if [ -n "${ALF_FLAGS_EXT+x}" ]; then
  printf "\nAppending additional compiler flag '%s'\n" "${ALF_FLAGS_EXT}"
fi

ALF_FLAGS_QRREF="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
# Modules need to know the programm configuration since entanglement needs MPI
ALF_FLAGS_MODULES="${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_FLAGS_EXT}"
ALF_FLAGS_ANA="${F90USEFULFLAGS} ${F90OPTFLAGS} ${ALF_INC} ${ALF_FLAGS_EXT}"
ALF_FLAGS_PROG="${F90USEFULFLAGS} ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_INC} ${ALF_FLAGS_EXT}"
# Control with flags -DHDF5 -DHDF5_ZLIB -DOBS_LEGACY, which observable format to use
if [ "${HDF5_ENABLED}" = "1" ]; then
  ALF_FLAGS_MODULES="${ALF_FLAGS_MODULES} ${INC_HDF5} -DHDF5 -DHDF5_ZLIB"
  ALF_FLAGS_ANA="${ALF_FLAGS_ANA} ${INC_HDF5} -DHDF5 -DHDF5_ZLIB"
  ALF_FLAGS_PROG="${ALF_FLAGS_PROG} -DHDF5 -DHDF5_ZLIB"
fi
export ALF_FLAGS_QRREF
export ALF_FLAGS_MODULES
export ALF_FLAGS_ANA
export ALF_FLAGS_PROG

rm -r "$tmpdir"
printf "\n${GREEN}Temporary directory %s deleted${NC}\n" "$tmpdir"

printf "\nTo compile your program use:    'make'\n\n"
