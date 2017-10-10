obj: bicgstab.c++ bicgstab.hpp ../dataparallel.d/dataparallel.c++
	rm -f lib/libbicgstab.a
	c++ -g -c bicgstab.c++ ../dataparallel.d/dataparallel.c++
	c++ -g -c ../TDMA.d/TDMA.c++
	ar cr lib/libbicgstab.a bicgstab.o dataparallel.o TDMA.o
	rm -f *.o

clean: ;
	rm -f *.o a.out

a.out: ;
	c++ sample.c++ -lbicgstab -lOpenCL
