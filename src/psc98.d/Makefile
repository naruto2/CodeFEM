cl_a.out: cl_main.c++ lib/libbicgstab.a
	cc -c yosen.c
	c++ cl_main.c++ cl_psc98.c++ yosen.o -lbicgstab -lOpenCL -o cl_a.out

clean: ;
	rm -f cl_a.out *.o lib/libpsc98.a

obj: ;
	cc -c yosen.c
	cc -c honsen.c
	c++ -c psc98.c++
	ar cr lib/libpsc98.a yosen.o honsen.o psc98.o
	rm *.o

time: ;
	export PSC98=0; time ./cl_a.out 120
	export PSC98=1; time ./cl_a.out 120
	export PSC98=2; time ./cl_a.out 120
	export PSC98=3; time ./cl_a.out 120
	export PSC98=4; time ./cl_a.out 120
	export PSC98=5; time ./cl_a.out 120
