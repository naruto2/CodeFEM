PROG =	glirulus
SRCS =	glirulus.c++ \
	op.c++ \
	mmio.c++ \
	GLU1_cpp.c++ 
CSRCS=	GLU1.c

OBJS =	$(SRCS:.c++=.o)
COBJS=	$(CSRCS:.c=.o)
.SUFFIXES: .c++ .c a.o

all: $(PROG)

$(PROG) : $(OBJS) $(COBJS)
	c++ $(OBJS) $(COBJS) -o $@ -lf2c

.c++.o:
	c++  -c $<

.c.o:
	cc -c $<


clean: ;
	rm -f glirulus *.o lib/libglirulus.a

obj: $(OBJS) $(COBJS)
	ar cr lib/libglirulus.a $(OBJS) $(COBJS)
