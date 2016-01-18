a.out: main.c++ solver.hpp obj
	cc -c yosen.c
	c++ -g main.c++ yosen.o -std=c++11 lib/libmatrix.a lib/libsolver.a

clean: ;
	rm -f a.out *.o lib/libsolver.a

obj: ;
	c++ -c blockmatrix.c++ -std=c++11
	ar cr lib/libsolver.a blockmatrix.o
