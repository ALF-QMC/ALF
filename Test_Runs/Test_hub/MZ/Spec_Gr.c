export ANNAL=${DIR}"/Analysis_8/"

{
declare -i nf
declare -i nfc

read nf 
let nfc=1
rm A_om.dat
while [ $nfc -le  $nf ];
    do
       read file
       export file1=Green_$file
       export file2=g_$file
       echo $file1
       cd $file1
       wc -l $file2 | sed s/${file}// > g_dat
       cat  $file2  >> g_dat       
       cp ../paramSAC .
       $ANNAL/Max_SAC.out
       while read line           
       do           
         echo $nfc" "$line   >> ../Green_om.dat
       done < Aom_ps_10
       echo  >> ../A_om.dat
       cd  ..
       let nfc=nfc+1
    done
rm ineq
} < Dyn_4

