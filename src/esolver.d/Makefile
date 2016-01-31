a.out: main.c++  obj
	c++  main.c++ -std=c++11 lib/*.a

clean: ;
	rm -f a.out *.o lib/libesolver.a

obj: ;
	c++ -c -std=c++11 minesolver.c++ maxesolver.c++
	ar cr lib/libesolver.a *.o
	rm *.o
