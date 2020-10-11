Name=(R_Om_st R_Om_en )

for File_sr in  "Green"  "SpinZ" "Den" ; do

case $File_sr in
     "Green")
     Variable=( -8.0 8.0 )
     Kpoints=( "0.00_0.00" "0.14_0.00" "0.27_0.00" "0.41_0.00" "0.55_0.00" "0.68_0.00" "0.82_0.00" "0.96_0.00" "1.09_0.00" \
         "1.23_0.00" "1.37_0.00" "1.50_0.00" "1.64_0.00" "1.78_0.00" "1.91_0.00" "2.05_0.00" "2.19_0.00" "2.32_0.00" \
	 "2.46_0.00" "2.60_0.00" "2.73_0.00" "2.87_0.00" "3.01_0.00" "3.14_0.00" )
     ;;
     "SpinZ")
     Variable=( 0.0 8.0 )
     Kpoints=( "0.14_0.00" "0.27_0.00" "0.41_0.00" "0.55_0.00" "0.68_0.00" "0.82_0.00" "0.96_0.00" "1.09_0.00" \
         "1.23_0.00" "1.37_0.00" "1.50_0.00" "1.64_0.00" "1.78_0.00" "1.91_0.00" "2.05_0.00" "2.19_0.00" "2.32_0.00" \
	 "2.46_0.00" "2.60_0.00" "2.73_0.00" "2.87_0.00" "3.01_0.00" "3.14_0.00" )
     ;;
     "Den")
     Variable=( 0.0 8.0 )
     Kpoints=( "0.00_0.00" "0.14_0.00" "0.27_0.00" "0.41_0.00" "0.55_0.00" "0.68_0.00" "0.82_0.00" "0.96_0.00" "1.09_0.00" \
         "1.23_0.00" "1.37_0.00" "1.50_0.00" "1.64_0.00" "1.78_0.00" "1.91_0.00" "2.05_0.00" "2.19_0.00" "2.32_0.00" \
	 "2.46_0.00" "2.60_0.00" "2.73_0.00" "2.87_0.00" "3.01_0.00" "3.14_0.00" )
esac	
rm  $File_sr"_Spectral.dat"
let n=0
for K in "${Kpoints[@]}"; do 
    export file=$File_sr"_"$K
    cd $file
    pwd
    cp ../parameters .
    let i=0
    while [  $i -lt 2 ]; do
	sed s/${Name[$i]}/${Variable[$i]}/  parameters  > tmp
	mv tmp parameters
	let i=i+1 
    done
    if [ ! -e "Green" ]; then
       ${ALF_DIR}/Analysis/Max_SAC.out
    fi
    cat Green | sed s/X/${n}/ >> ../${File_sr}"_Spectral.dat"
    cd ..
    echo >> $File_sr"_Spectral.dat"
    echo >> $File_sr"_Spectral.dat"
    let n=n+1
done

export Plot_File=${File_sr}"_Spectral.gnu"
echo $Plot_File
if [ ! -e ${Plot_File} ]; then
   cp Plot_Spectral_template.gnu  ${Plot_File}
   let i=1
   while [  $i -lt 3 ]; do
     sed s/${Name[$i]}/${Variable[$i]}/  $Plot_File  > tmp
     mv tmp $Plot_File
     let i=i+1 
   done	
   sed s/"R_Plot_Name"/${File_sr}/  $Plot_File  > tmp
   mv tmp $Plot_File
   let n=n-1
   sed s/"R_xrange"/$n/  $Plot_File  > tmp
   mv tmp $Plot_File
   sed s/"R_Plot_File"/$File_sr"_Spectral.dat"/  $Plot_File  > tmp
   mv tmp $Plot_File
fi

done
