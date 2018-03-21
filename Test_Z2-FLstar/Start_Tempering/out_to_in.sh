declare -i n
declare -i n1

n1=`ls confout_* | wc -l`

let n=0
let n1=n1-1
while [ $n -le  $n1 ];
   do
     file_unzipped="confout_"$n
     export file_out="confout_"$n".gz"
     if [ -f "$file_unzipped" ]
     then
       gzip -9 "$file_unzipped"
     fi
     export file_in="confin_"$n".gz"
     echo $file_out  $file_in
     mv $file_out $file_in
     let n=n+1
   done
