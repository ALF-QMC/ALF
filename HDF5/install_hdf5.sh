#!/usr/bin/env bash
# Script for automatically downloading and installing HDF5 in current directory
# Needs the following environment variables:
#   CC: C compiler
#   FC: Fortran compiler
#   CXX: C++ compiler
#   HDF5_DIR: Diretory, in which HDF5 gets installed
H5_major="$1"
H5_minor="$2"
H5_patch="$3"
H5_suff="$4"

if [ -d "$HDF5_DIR" ]; then
  printf "\e[31mDirectory %s already exists, aborting HDF5 installation.\e[0m\n" "$HDF5_DIR" 1>&2
  exit 1
fi

command -v curl > /dev/null
CURL_AVAIL=$?
command -v wget > /dev/null
WGET_AVAIL=$?

if [ $CURL_AVAIL -ne 0 ] && [ $WGET_AVAIL -ne 0 ]; then
  printf "\e[31m==== Neither curl nor wget available!                   =====\e[0m\n" 1>&2
  printf "\e[31m==== One of the two is required to download HDF5 source =====\e[0m\n" 1>&2
  exit 1
fi

if ! command -v "$CC" > /dev/null; then
  printf "\e[31m==== C compiler <%s> not available =====\e[0m\n" "$CC" 1>&2
  exit 1
fi
if ! command -v "$FC" > /dev/null; then
  printf "\e[31m==== FORTRAN compiler <%s> not available =====\e[0m\n" "$FC" 1>&2
  exit 1
fi
if ! command -v "$CXX" > /dev/null; then
  printf "\e[31m==== C++ compiler <%s> not available =====\e[0m\n" "$CXX" 1>&2
  exit 1
fi

# Create temporary directory
tmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir')
printf "\033[0;32mTemporary directory %s created\e[0m\n" "$tmpdir"
cd "$tmpdir" || exit 1

printf "\033[0;32m========== Downloading HDF5 source ==========\e[0m\n" 1>&2

H5_SRC="https://support.hdfgroup.org/releases/hdf5/v${H5_major}_${H5_minor}/v${H5_major}_${H5_minor}_${H5_patch}/downloads/hdf5-${H5_major}.${H5_minor}.${H5_patch}${H5_suff}.tar.gz"
echo "From $H5_SRC"
source_dir="hdf5-${H5_major}.${H5_minor}.${H5_patch}${H5_suff}"
if [ $CURL_AVAIL -eq 0 ]; then
  curl "$H5_SRC" | tar xz || exit 1
else
  wget -O- "$H5_SRC" | tar xz || exit 1
fi

export CC FC CXX
printf "\033[0;32m=== Build with the following compilers C: %s, Fortran: %s, C++: %s \e[0m\n" "$CC" "$FC" "$CXX" 1>&2

"$source_dir/configure" --prefix="$HDF5_DIR" --libdir="$HDF5_DIR/lib" --enable-fortran --enable-shared=no --enable-tests=no
if ! make; then
  printf "\e[31m=== Compilation with compilers %s %s in directory %s failed ===\e[0m\n" "$CC" "$FC" "$PWD" 1>&2
  rm -rf "$HDF5_DIR"
  exit 1
fi
#make check
if ! make install; then
  printf "\e[31m=== Installation of HDF5 in directory %s failed ===\e[0m\n" "$HDF5_DIR" 1>&2
  rm -rf "$HDF5_DIR"
  exit 1
fi
#make check-install

printf "\033[0;32mYou can delete the temporary directory %s\e[0m\n" "$tmpdir"
