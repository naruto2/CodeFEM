a.out: ;
	cc *.c -lf2c -lm
	./a.out > cout.txt
	diff cout.txt fout.txt
	rm -f a.out

fout: ;
	f77 *.f
	./a.out > fout.txt
	rm -f a.out

clean: ;
	rm -f a.out

psc: ;
	cc -c PCGS.c
	c++ -c PCGS_cpp.c++ psc98.c++
	c++ psc98.o PCGS_cpp.o PCGS.o lib/*.a -lf2c -lm
