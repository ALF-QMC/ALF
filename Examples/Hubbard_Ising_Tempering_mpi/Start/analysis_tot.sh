for filename in Temp_*; do
    cd $filename
    cp ../analysis.sh .
    bash analysis.sh
    cd ..
done

