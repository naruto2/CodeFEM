PROG =  a.out
SRCS =	ncpolynomial1.c++ \
	printmatrix.c++ \
	plotncpolynomial1.c++ \
	squaremesh.c++ \
	batonth.c++ \
	cont.c++ \
	stokes.c++


OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -std=gnu++11
LDFLAGS = lib/*.a -lOpenCL

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS) main.c++
	c++ main.c++ -std=gnu++11 $(OBJS) -o $@ $(LDFLAGS) -pthread

.c++.o:
	c++ $(CPPFLAGS) -c $< 

obj: $(OBJS)
	ar cr lib/libstokes.a $(OBJS)
	rm *.o

clean :
	rm -f $(OBJS) $(PROG) lib/libstokes.a

test : $(PROG)
	./a.out 10x10.xyc

world: ;
	(cd ../../bin; ./mklink)
