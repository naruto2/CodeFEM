a.out: main.c++
	c++ -Wall main.c++ -Llib -lpsc98 -lViennaCL -lop -lm  -lbicgstab -lTDMA -lOpenCL

amg: amg.c++
	c++ amg.c++ -Llib -lpsc98 -o amg

time: a.out
	export PSC98=0; time ./a.out
	export PSC98=1; time ./a.out
	export PSC98=2; time ./a.out
	export PSC98=3; time ./a.out
	export PSC98=4; time ./a.out
	export PSC98=5; time ./a.out

honsen: main.c++
	cc -c honsen.c
	c++ -Wall main.c++ honsen.o -lOpenCL -lm -o honsen

honsenrun: honsen
	export PSC98=1; time ./honsen
	export PSC98=2; time ./honsen
	export PSC98=3; time ./honsen
	export PSC98=4; time ./honsen
	export PSC98=5; time ./honsen

clean: ;
	rm -f a.out amg honsen amg.o yosen.o honsen.o ViennaCL.o lib/libViennaCL.a 

obj: ;
	c++ -c ViennaCL.c++
	ar cr lib/libViennaCL.a ViennaCL.o
	rm -f ViennaCL.o
