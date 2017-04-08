Variable=(Run L_R  B_R  J_R Jz_R  U_R Nsweeps_R  Nbins_R Nwrap_R N_jobs)
Name=(Run L_R  B_R  J_R Jz_R  U_R Nsweeps_R  Nbins_R Nwrap_R N_jobs)

{
read -a Variable
echo ${Variable[@]}
echo ${Variable[0]}
while [ ! ${Variable[0]} = "stop" ];   do
    if [ ${Variable[0]} = "Y" ]; then
     

    export B_R_dir=`echo ${Variable[2]} | sed s/"\.0"//`
    export J_R_dir=`echo ${Variable[3]} | sed s/"\.0"//`
    export Jz_R_dir=`echo ${Variable[4]}  | sed s/"\.0"//`
    export U_R_dir=`echo ${Variable[5]}  | sed s/"\.0"//`
    
    export Dir="L"${Variable[1]}"_U"$U_R_dir"_B"$B_R_dir"_J"$J_R_dir"_Jp"$Jz_R_dir
    echo $Dir
    if [ ! -e $Dir ]; then
      mkdir $Dir
      cd  $Dir
      cp  ../Start/* . 
      cd ..
    else
      cd $Dir
      if [ -e  "confout_0" ]; then
          cp ../Start/out_to_in.c .
          bash out_to_in.c
      fi
      cp ../Start/parameters .
      cd ..
    fi
    cd $Dir

    i=1
    while [  $i -lt 9 ]; do
	sed s/${Name[$i]}/${Variable[$i]}/    parameters  > tmp
	mv tmp parameters
        let i=i+1 
    done

    #sed s/Dir_r/$Dir/ job.q  > job1.q
    #mv job1.q job.q

    #export Dir1="K"${Variable[1]}"_h"$H_R_dir
    #export Run_name=$Dir1
    #echo $Run_name
    #sed s/Run_name/$Run_name/ chain_jobs.c > job1.q
    #sed s/Nj_r/${Variable[9]}/         job1.q > chain_jobs.c
    #if [ ${Variable[8]} -gt 0 ]; then 
    #	bash chain_jobs.c
    #fi

   cd ..
   fi
read -a Variable
done
} <  Sims_Jp01
