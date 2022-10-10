#
for th in 64 32 16 8 4 2 1; do
	export ACC_NUM_CORES=$th
	echo "Running threads = "$ACC_NUM_CORES
	./bgk_thermal > out.$ACC_NUM_CORES.log
	md5sum Output*50000.vtk > check.$ACC_NUM_CORES.log
done
