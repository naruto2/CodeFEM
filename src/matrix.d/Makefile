a.out: matrix.hpp main.c++ matrix.c++ 
	c++ -Wall main.c++ matrix.c++ -std=c++11

clean: ;
	rm -f a.out lib/libmatrix.a

obj: ;
	c++ -c  matrix.c++ -std=c++11
	ar cr lib/libmatrix.a *.o
	rm *.o
