obj: bicgstab.c++ bicgstab.hpp ../dataparallel.d/dataparallel.c++
	rm -f lib/libbicgstab.a
	c++ -g -c bicgstab.c++ ../dataparallel.d/dataparallel.c++
	ar cr lib/libbicgstab.a bicgstab.o dataparallel.o
	rm -f *.o

clean: ;
	rm -f *.o a.out

a.out: ;
	c++ sample.c++ -lbicgstab -lOpenCL
