a.out: main.c++ yosen.o 
	c++ main.c++ yosen.o -lOpenCL

yosen.o: yosen.c
	cc -c yosen.c

clean: ;
	rm -f a.out *.o lib/libEigen.a

obj:  Elu.c++
	c++ -c Elu.c++
	ar cr lib/libEigen.a Elu.o
	rm -f Elu.o

