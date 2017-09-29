a.out: main.c++ 
	c++ main.c++ lib/libstwart.a -lf2c

clean: ;
	rm -f a.out lib/libstwart.a

obj: FORFUL.c MARKWZ.c MARKWZT.c SINGLE.c STWART.c stwart_cpp.c++
	cc -c *.c
	c++ -c stwart_cpp.c++
	ar cr lib/libstwart.a *.o
	rm *.o

