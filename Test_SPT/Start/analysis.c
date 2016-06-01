export ANNAL=${DIR}"/Analysis_7/"

$ANNAL/jackv5.out 
for filename in *_eq; do
    echo $filename
    export filename1=$filename"J"
    if [ "$filename1" -ot "$filename" ]; then
       cp $filename ineq
       $ANNAL/cov_eq.out
       mv "equalJ"  $filename"J"
    fi
done


#for filename in *_tau; do
#    echo $filename
#    cp $filename intau
#    $ANNAL/cov_tau.out
#    for filename1 in g_*; do
#	echo $filename1
#	export Name=`echo ${filename} | sed s/_tau//`	
#	export Dir=`echo ${filename1} | sed s/g_/${Name}_/`
#	echo $Dir
#	mkdir $Dir
#	cd  $Dir 
#	mv ../$filename1 .
#	cd  ..
#    done
#
#done

