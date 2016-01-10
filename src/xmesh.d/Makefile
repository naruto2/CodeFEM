PROG =  xmesh
SRCS =	rotate_right.c++ \
	rotate_left.c++ \
	N_set.c++ \
	degeneracy.c++ \
	super_triangle.c++ \
	LawsonSwap.c++ \
	vanish_boundary_triangle.c++ \
	delaunay.c++ \
	count_nodes.c++ \
	generate_fN.c++ \
	outmesh.c++ \
	plotmesh.c++ \
	in2xyc.c++ 


OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -Wall
LDFLAGS = lib/*.a

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS)
	c++ main.c++ $(OBJS) -o $@ $(LDFLAGS)
	rm $(OBJS)

.c++.o:
	c++ $(CPPFLAGS) -c $< 

obj: $(OBJS)
	ar cr lib/libxyc2msh.a $(OBJS)
	rm -f $(OBJS)

clean :
	rm -f $(OBJS) $(PROG) lib/libxyc2msh.a kanto.msh1 pt.msh1 test

test : $(PROG)
	./xmesh kanto.xyc -o kanto.msh1
	./xmesh pt.xyc -o pt.msh1
	./xmesh .msh0 kanto.msh1
	diff pt.msh0 pt.msh1
	touch test
