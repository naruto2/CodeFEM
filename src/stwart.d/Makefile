a.out: main.c++ lib/libstwart.a
	c++ main.c++ lib/libstwart.a -lf2c
	./a.out

obj: ;
	cc -c *.c
	c++ -c stwart_cpp.c++
	ar cr lib/libstwart.a *.o
	rm *.o

clean: a.out
	rm a.out lib/libstwart.a
