a.out: main.c++
	(cd ../bicgstab.d; make obj)
	c++ -g main.c++ -lnavierstokes -lbicgstab -lxmesh -lforeach -lOpenCL

clean: ;
	rm -f a.out lib/libnavierstokes.a

obj: ;
	c++ -c navierstokes.c++ 
	ar cr lib/libnavierstokes.a navierstokes.o
	rm -f navierstokes.o
