#
for inp in 0 1 2 3 ; do
	echo "Running size = input.dat."$inp
	cp input.dat.$inp input.dat
	./bgk_thermal > out.size.$inp.log
done
