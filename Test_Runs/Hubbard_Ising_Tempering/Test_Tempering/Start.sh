Variable=(H_R  )
Name=(H_R )

{
read -a Variable
let n=0
while [ ! ${Variable[0]} = "stop" ];   do
     mkdir Temp_$n
     cd  Temp_$n
     cp ../parameters . 
     sed s/${Name[0]}/${Variable[0]}/  parameters  > tmp
     mv tmp parameters
     cp ../seeds .
     cp  ../out_to_in.sh .
     cd ..
     let n=n+1
read -a Variable
done
} < Sims

