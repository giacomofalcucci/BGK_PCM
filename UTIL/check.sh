for time in 00000 05000 10000 15000 20000 30000 40000 50000; do 
    md5sum STEP_CUF/REFERENCE_1/Output_000$time.vtk 
    md5sum STEP_CUF/RUN/Output_000$time.vtk 
    echo "t = ", $time, " (press enter to continue)"
    read
done

