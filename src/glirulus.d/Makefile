PROG =	glirulus
SRCS =	glirulus.c++ \
	mmio.c++ \
	GLU1_cpp.c++ 
CSRCS=	GLU1.c

OBJS =	$(SRCS:.c++=.o)
COBJS=	$(CSRCS:.c=.o)
.SUFFIXES: .c++ .c a.o

all: $(PROG)

$(PROG) : $(OBJS) $(COBJS)
	c++ main.c++ $(OBJS) $(COBJS) -o $@ -Llib -lop -lf2c -lbicgstab -lOpenCL -lEigen -lViennaCL -lTDMA -lperfectpivot -ljacobi
	rm *.o

.c++.o:
	c++  -c $<

.c.o:
	cc -c $<


clean: ;
	rm -f glirulus *.o lib/libglirulus.a

obj: $(OBJS) $(COBJS)
	ar cr lib/libglirulus.a $(OBJS) $(COBJS)
	rm *.o
