a.out: main.c++ 
	c++ main.c++ lib/libmatrix.a lib/libsolver.a lib/libpsc98.a

clean: ;
	rm -f a.out *.o lib/libpsc98.a

obj: ;
	cc -c yosen.c
	cc -c honsen.c
	c++ -c psc98.c++
	ar cr lib/libpsc98.a yosen.o honsen.o psc98.o

