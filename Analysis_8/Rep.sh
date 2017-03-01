for filename in *.f90; do
  echo $filename
  mv $filename file.f90
  sed s/Kind=8/"Kind=Kind(0.d0)"/ file.f90 > file1.f90
  sed s/kind=8/"Kind=Kind(0.d0)"/ file1.f90 > file.f90
  sed s/KIND=8/"Kind=Kind(0.d0)"/ file.f90 > file1.f90
  sed s/"=double"/"=Kind(0.d0)"/ file1.f90 > file.f90
  mv file.f90 $filename
done

