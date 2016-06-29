PROG =  a.out
SRCS =	Rf.c++ \
	ary.c++ \
	foreach.c++ \
	ncpolynomial1.c++ \
	forgammap1.c++ \
	printmatrix.c++ \
	plotncpolynomial1.c++ \
	squaremesh.c++

OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -std=gnu++11 -Wall -fpermissive -I. 
LDFLAGS = lib/*.a -lOpenCL

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS) main.c++
	c++ main.c++ -std=c++11 $(OBJS) -o $@ $(LDFLAGS) 

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
