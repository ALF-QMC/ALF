declare -i n
declare -i n1

n1=`ls confout_* | wc -l`

let n=0
let n1=n1-1
while [ $n -le  $n1 ];
   do
     export file_out="confout_"$n
     export file_in="confin_"$n
     echo $file_out  $file_in
     mv $file_out $file_in
     let n=n+1
   done
