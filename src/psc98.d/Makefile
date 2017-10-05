cl_a.out: main.c++ lib/libbicgstab.a
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
