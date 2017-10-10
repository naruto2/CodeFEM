cl_a.out: cl_main.c++ lib/libbicgstab.a lib/libpsc98.a
	(cd ../bicgstab.d; make obj)
	c++ cl_main.c++  -lpsc98 -lbicgstab -lOpenCL -o cl_a.out

clean: ;
	rm -f cl_a.out *.o lib/libpsc98.a

lib/libpsc98.a: obj
	echo -n

obj: ;
	cc -c yosen.c
	c++ -c psc98.c++
	c++ -c cl_psc98.c++
	ar cr lib/libpsc98.a yosen.o psc98.o cl_psc98.o
	rm *.o

time: ;
	export PSC98=0; time ./cl_a.out 150
	export PSC98=1; time ./cl_a.out 150
	export PSC98=2; time ./cl_a.out 150
	export PSC98=3; time ./cl_a.out 150
	export PSC98=4; time ./cl_a.out 150
	export PSC98=5; time ./cl_a.out 150

lib/libbicgstab.a: ;
	(cd ../bicgstab.d; make obj)
