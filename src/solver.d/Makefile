a.out: main.c++ solver.hpp obj
	c++ -g main.c++ -std=c++11 lib/libmatrix.a lib/libsolver.a lib/libpsc98.a

clean: ;
	rm -f a.out *.o lib/libsolver.a

obj: ;
	c++ -c blockmatrix.c++ -std=c++11
	c++ -c solver.c++ -std=c++11
	ar cr lib/libsolver.a blockmatrix.o solver.o
	rm *.o
