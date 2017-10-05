a.out: main.c++ lib/libbicgstab.a
	c++ main.c++  lib/libbicgstab.a lib/libxmesh.a lib/libforeach.a -lOpenCL


lib/libbicgstab.a: ;
	(cd ../bicgstab.d/;make)

clean: ;
	rm -f a.out lib/libbicgstab.a

obj: ;
	echo -n
