a.out: main.c++ 
	c++ -g maxesolver.c++ main.c++ -std=c++11 lib/libmatrix.a lib/libsolver.a lib/libpsc98.a

clean: ;
	rm -f a.out *.o 

obj: ;
	c++ -c blockmatrix.c++ -std=c++11
	c++ -c solver.c++ -std=c++11
	ar cr lib/libsolver.a blockmatrix.o solver.o
