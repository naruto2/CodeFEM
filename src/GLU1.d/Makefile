a.out: main.c++ 
	c++ -g main.c++ lib/libGLU1.a -lf2c

clean: ;
	rm -f a.out lib/libGLU1.a *.o

obj: GLU1.c GLU1_cpp.c++
	cc -g -c *.c
	c++ -g -c *.c++
	ar cr lib/libGLU1.a *.o
	rm *.o

