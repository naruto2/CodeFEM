a.out: main.c++ dataparallel.c++ cl_bicgstab_kernel.cl
	c++ main.c++ dataparallel.c++ -lOpenCL lib/libop.a
	./a.out

clean: ;
	rm -f ./a.out

obj: ;
	echo -n
