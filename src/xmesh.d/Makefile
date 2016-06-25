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
	in2xyc.c++ \
	sortmesh.c++


OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -Wall
LDFLAGS = lib/*.a

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS)
	c++ main.c++ $(OBJS) -o $@ $(LDFLAGS)


.c++.o:
	c++ $(CPPFLAGS) -c $< 

obj: $(OBJS)
	ar cr lib/libxmesh.a $(OBJS)
	rm *.o

clean :
	rm -f $(OBJS) $(PROG) lib/libxmesh.a kanto.msh pt.msh test

test : $(PROG)
	rm -f kanto.msh pt.msh
	./xmesh kanto.xyc -o kanto.msh
	./xmesh pt.xyc -o pt.msh
	diff kanto_org.msh kanto.msh
	diff pt_org.msh pt.msh
	touch test
