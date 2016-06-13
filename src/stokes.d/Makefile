PROG =  a.out
SRCS =	Rf.c++ \
	ary.c++ \
	foreach.c++ \
	ncpolynomial1.c++

OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -Wall -fpermissive
LDFLAGS = lib/*.a

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS) main.c++
	c++ main.c++ $(OBJS) -o $@ $(LDFLAGS) 

.c++.o:
	c++ $(CPPFLAGS) -c $< 

obj: $(OBJS)
	ar cr lib/libstokes.a $(OBJS)
	rm *.o

clean :
	rm -f $(OBJS) $(PROG) 

test : $(PROG)
	./a.out 10x10.xyc
