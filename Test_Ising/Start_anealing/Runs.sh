Variable=(Run h_R )
Name=(Run h_R  )

{
read -a Variable
echo ${Variable[@]}
echo ${Variable[0]}
while [ ! ${Variable[0]} = "stop" ];   do
    if [ ${Variable[0]} = "Y" ]; then

    cp parameters_hvar parameters

    i=1
    while [  $i -lt 2 ]; do
        sed s/${Name[$i]}/${Variable[$i]}/    parameters  > tmp
        mv tmp parameters
        let i=i+1 
    done

    ../../Prog_8/Hubbard_Ising.out
    bash analysis.sh
    export tmp=`grep "OBS :    1 " Flux_scalJ | sed s/"OBS :    1"//`
    echo ${Variable[1]}"  "$tmp >> Flux_sweepsJ
    rm *_tau *_scal *_eq
    bash out_to_in.sh
    fi
read -a Variable
done
} < Sims
