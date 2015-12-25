a.out: main.c++
	cc -c yosen.c
	c++ main.c++ yosen.o -lOpenCL -lm

run: a.out
	export PSC98=1; time ./a.out
	export PSC98=2; time ./a.out
	export PSC98=3; time ./a.out
	export PSC98=4; time ./a.out
	export PSC98=5; time ./a.out

honsen: main.c++
	cc -c honsen.c
	c++ main.c++ honsen.o -lOpenCL -lm -o honsen

honsenrun: honsen
	export PSC98=1; time ./honsen
	export PSC98=2; time ./honsen
	export PSC98=3; time ./honsen
	export PSC98=4; time ./honsen
	export PSC98=5; time ./honsen

clean: ;
	rm -f a.out honsen yosen.o honsen.o

obj: ;
	touch lib
