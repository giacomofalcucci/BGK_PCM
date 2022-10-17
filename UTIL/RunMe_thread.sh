#
for th in 64 32 16 8 4 2 1; do
	export ACC_NUM_CORES=$th
	export OMP_NUM_THREADS=$th
	echo "Running threads = "$ACC_NUM_CORES, $OMP_NUM_THREADS
	./bgk_thermal > out.$th.log
	md5sum Output*50000.vtk > check.$th.log
done
