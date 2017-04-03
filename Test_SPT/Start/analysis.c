export ANNAL=${HOME}"/Documents/Doktorarbeit/General_QMCT_code/Analysis_8"

$ANNAL/jackv5.out 
for filename in *_eq; do
    echo $filename
    export filename1=$filename"J"
    if [ "$filename1" -ot "$filename" ]; then
       cp $filename ineq
       $ANNAL/cov_eq.out
       rm ineq
       mv "equalJ"  $filename"J"
       mv "equalJR"  $filename"RJ"
    fi
done

for filename in *_sus; do
    echo $filename
    export filename1=$filename"J"
    if [ "$filename1" -ot "$filename" ]; then
       cp $filename ineq
       $ANNAL/cov_eq.out
       rm ineq
       mv "equalJ"  $filename"J"
       mv "equalJR"  $filename"RJ"
    fi
done


for filename in *_tau; do
    echo $filename
    cp $filename intau
    $ANNAL/cov_tau.out
    rm intau
    filetmp=$filename"_timeline"
    mv timeline $filetmp
    for filename1 in g_*; do
	echo $filename1
	export Name=`echo ${filename} | sed s/_tau//`	
	export Dir=`echo ${filename1} | sed s/g_/${Name}_/`
	echo $Dir
	mkdir -p $Dir
	cd  $Dir 
	mv ../$filename1 .
	cd  ..
    done

done

