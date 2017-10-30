psc98: psc98.c++ lib/libbicgstab.a lib/libpsc98.a
	c++ psc98.c++  -Llib -lpsc98 -lOpenCL -lglirulus -lViennaCL -lbicgstab -lEigen -lTDMA -lop -lperfectpivot -o psc98

clean: ;
	rm -f psc98 *.o lib/libpsc98.a

lib/libpsc98.a: obj
	echo -n

obj: ;
	cc -c yosen.c
	c++ -c cl_psc98.c++
	ar cr lib/libpsc98.a yosen.o cl_psc98.o
	rm *.o

time: ;
	export PSC98=0; time ./psc98 -solver vcl_gmres
	export PSC98=1; time ./psc98 -solver vcl_gmres
	export PSC98=2; time ./psc98 -solver vcl_gmres
	export PSC98=3; time ./psc98 -solver vcl_gmres
	export PSC98=4; time ./psc98 -solver vcl_gmres
	export PSC98=5; time ./psc98 -solver vcl_gmres

lib/libbicgstab.a: ;
	(cd ../bicgstab.d; make obj)
