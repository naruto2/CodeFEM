a.out: main.c++
	c++ main.c++ -lnavierstokes -lbicgstab -lxmesh -lforeach -lOpenCL

clean: ;
	rm -f a.out

obj: ;
	c++ -c navierstokes.c++ 
	ar cr lib/libnavierstokes.a navierstokes.o
	rm -f navierstokes.o
