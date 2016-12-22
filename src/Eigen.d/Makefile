a.out: main.c++ yosen.o 
	c++ main.c++ yosen.o -lOpenCL

yosen.o: yosen.c
	cc -c yosen.c

clean: ;
	rm -f a.out *.o
