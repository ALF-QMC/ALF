let n=0
for filename in Temp_*; do
    cd $filename
    bash out_to_in.sh
    let n=n+1
    cd ..
done

mpirun -np $n ${BSS}/Prog_8/SPT.out
