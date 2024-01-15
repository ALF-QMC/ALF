#!/bin/bash
shopt -s nocasematch 

Name=(R_Om_st R_Om_en R_Ndis R_Stochastic R_N_Cov)
export New_run="True"

#for File_sr in  "Green" "Den" "SpinZ"  ; do
for File_sr in  "Green" ; do

case $File_sr in
     "Green")
     Variable=( -8.0 8.0 600 False 1)
     Kpoints=( "0.00_0.00" "0.14_0.00" "0.27_0.00" "0.41_0.00" "0.55_0.00" "0.68_0.00" "0.82_0.00" "0.96_0.00" "1.09_0.00" \
         "1.23_0.00" "1.37_0.00" "1.50_0.00" "1.64_0.00" "1.78_0.00" "1.91_0.00" "2.05_0.00" "2.19_0.00" "2.32_0.00" \
	     "2.46_0.00" "2.60_0.00" "2.73_0.00" "2.87_0.00" "3.01_0.00" "3.14_0.00" )
     ;;
     "SpinZ")
     Variable=( 0.0 8.0 600 False 1)
     Kpoints=( "0.14_0.00" "0.27_0.00" "0.41_0.00" "0.55_0.00" "0.68_0.00" "0.82_0.00" "0.96_0.00" "1.09_0.00" \
         "1.23_0.00" "1.37_0.00" "1.50_0.00" "1.64_0.00" "1.78_0.00" "1.91_0.00" "2.05_0.00" "2.19_0.00" "2.32_0.00" \
	     "2.46_0.00" "2.60_0.00" "2.73_0.00" "2.87_0.00" "3.01_0.00" "3.14_0.00" )
     ;;
     "Den")
     Variable=( 0.0 8.0 600 False 1)
     Kpoints=(  "0.14_0.00" "0.27_0.00" "0.41_0.00" "0.55_0.00" "0.68_0.00" "0.82_0.00" "0.96_0.00" "1.09_0.00" \
         "1.23_0.00" "1.37_0.00" "1.50_0.00" "1.64_0.00" "1.78_0.00" "1.91_0.00" "2.05_0.00" "2.19_0.00" "2.32_0.00" \
	     "2.46_0.00" "2.60_0.00" "2.73_0.00" "2.87_0.00" "3.01_0.00" "3.14_0.00" )
esac	

if [[ "${Variable[3]}" == "true" ]] ; then
     export File_end="_Spectral.dat"
     export File_end_plot="_Spectral.gnu"
     export Green_file="Green"
     export Plot_name="$File_sr"_Spectral
     export Plot_Title="$File_sr"
else
     export File_end="_Spectral_cl.dat"
     export File_end_plot="_Spectral_cl.gnu"
     export Green_file="Green_cl"
     export Plot_name="$File_sr"_Spectral_cl
     export Plot_Title="$File_sr"_cl
fi
export Data_file=$File_sr$File_end
if [ -f $Data_file ]; then
    echo "Data File exists ${Data_file},  will  remove it"
    rm $Data_file
fi
(( n=0 ))
for K in "${Kpoints[@]}"; do 
    export file=$File_sr"_"$K
    cd "$file" || exit 
    cp ../parameters .
    (( i=0 ))
    while [ $i -lt 5 ]; do
	    sed s/"${Name[$i]}"/"${Variable[$i]}"/  parameters  > tmp
	    mv tmp parameters
	    (( i=i+1 )) 
    done
    if [[ "${New_run}" == "true" ]]; then
        rm $Green_file
    fi
    if [ ! -e "$Green_file" ]; then
        #echo "Run_code_here"
        "${ALF_DIR}"/Analysis/Max_SAC.out
    fi
    sed s/X/${n}/ "$Green_file"  >> ../${Data_file} 
    cd ..
    echo >> $Data_file
    echo >> $Data_file
    (( n=n+1 ))
done

export Plot_File=${File_sr}${File_end_plot}
#echo $Plot_File
if [ ! -e ${Plot_File} ]; then
   cp Plot_Spectral_template.gnu  ${Plot_File}
   (( i=0 ))
   while [  $i -lt 2 ]; do
     sed s/"${Name[$i]}"/"${Variable[$i]}"/  $Plot_File  > tmp
     mv tmp $Plot_File
     (( i=i+1 ))
   done	
   sed s/"R_Plot_Name"/${Plot_name}/  $Plot_File  > tmp
   mv tmp $Plot_File
   (( n=n-1 ))
   sed s/"R_xrange"/$n/  $Plot_File  > tmp
   mv tmp $Plot_File
   sed s/"R_Plot_File"/$Data_file/  "$Plot_File"  > tmp
   mv tmp $Plot_File
   sed s/"R_Plot_Title"/$Plot_Title/  "$Plot_File"  > tmp
   mv tmp $Plot_File
fi

done
