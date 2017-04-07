for filename in Temp_*; do
    cd $filename
    bash out_to_in.sh
    cd ..
done

