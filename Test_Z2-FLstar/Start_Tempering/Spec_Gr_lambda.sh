export ANNAL=$BSS"/Analysis/"

{
rm -f Green_lambda_om.dat
while read lambda dump 
  do
    dir="l_"$lambda
    cd $dir
    file="3.14_3.14"
    export file1=Green_$file
    export file2=g_$file
    echo $dir/$file1
    cd $file1
    wc -l $file2 | sed s/${file}// > g_dat
    cat  $file2  >> g_dat       
    cp ../../../paramSAC .
    $ANNAL/Max_SAC.out
    while read line           
    do           
      echo "1 "$line   >> ../Green_om.dat
      echo $lambda" "$line   >> ../../Green_lambda_om.dat
    done < Aom_ps_10
    echo "" >> ../../Green_lambda_om.dat
    cd  ../..
  done
rm -f ineq
} < DerivativeOfF
