a.out: op.hpp op.c++ main.c++
	c++ -Wall *.c++

obj: op.hpp op.c++
	c++ -Wall -c op.c++
	ar cr lib/libop.a *.o
	rm *.o


test: a.out
	./a.out -v aaa -f bbb

clean:
	rm -f a.out lib/libop.a *.o
