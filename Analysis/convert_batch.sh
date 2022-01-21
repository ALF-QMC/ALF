#!/bin/bash
# Script for converting plain text bins to HDF5 format.

function convert {
  echo "=== Converting $PWD ==="
  >&2 echo "=== Converting $PWD ==="
  find . -maxdepth 1 -name "*_scal" -exec "${dir_bin}/convert_scal.out" {} data.h5 \;
  find . -maxdepth 1 -name "*_eq"   -exec "${dir_bin}/convert_latt.out" {} data.h5 \;
  find . -maxdepth 1 -name "*_tau"  -exec "${dir_bin}/convert_latt.out" {} data.h5 \;
}
export -f convert

export dir_bin=$(dirname "$0")

# File to recognize in which folders to convert
ref="Ener_scal"

find . -name "$ref" -execdir bash -c 'convert' \;

