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
	sortmesh.c++ \
	p2.c++ \
	f2mesh.c++ \
	dimp2.c++ \
	delta.c++


OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS =
LDFLAGS = -Llib -lop

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS)
	c++ xmesh.c++ $(OBJS) -o $@ $(LDFLAGS)

.c++.o:
	c++ $(CPPFLAGS) -c $< 

obj: $(OBJS)
	ar cr lib/libxmesh.a $(OBJS)
	rm *.o

clean :
	rm -f $(OBJS) $(PROG) lib/libxmesh.a 


