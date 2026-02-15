#!/bin/bash
shopt -s extglob

out_to_in()
{
  if [ -e confout_0.h5 ]  && [ -e confout_0 ]; then
    >&2 echo "out_to_in.sh: Both plain text and HDF5 configuration present. Please remove one of them."
    exit 1
  fi
  
  if [ -e confout_0.h5 ]; then
    for file_out in confout_+([0-9]).h5; do
      file_in=${file_out/out/in}
      echo "mv $file_out to $file_in"
      mv "$file_out" "$file_in"
    done
  elif [ -e confout_0 ]; then
    for file_out in confout_+([0-9]); do
      file_in=${file_out/out/in}
      echo "mv $file_out to $file_in"
      mv "$file_out" "$file_in"
    done
  fi
}

for i in Temp_*; do
  (cd "$i" && out_to_in)
done
