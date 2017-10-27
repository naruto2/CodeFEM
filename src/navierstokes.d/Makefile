navierstokes: main.c++
	c++ -g main.c++ -Llib -lnavierstokes -lbicgstab -lxmesh -lforeach -lOpenCL -lop -lglirulus -lViennaCL -lEigen -lTDMA -lf2c -o navierstokes

clean: ;
	rm -f navierstokes lib/libnavierstokes.a

obj: ;
	c++ -c navierstokes.c++ 
	ar cr lib/libnavierstokes.a navierstokes.o
	rm -f navierstokes.o
