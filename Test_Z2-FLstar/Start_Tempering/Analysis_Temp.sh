for filename in Temp_*; do
    cd $filename
    bash ../analysis_Auto.sh
    cd ..
done
