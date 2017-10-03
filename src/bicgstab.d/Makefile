lib/libbicgstab.a: bicgstab.c++ bicgstab.hpp ../dataparallel.d/dataparallel.c++
	rm -f lib/libbicgstab.a
	c++ -c bicgstab.c++ dataparallel.c++
	ar cr lib/libbicgstab.a bicgstab.o dataparallel.o
