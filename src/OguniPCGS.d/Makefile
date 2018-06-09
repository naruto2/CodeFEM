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
