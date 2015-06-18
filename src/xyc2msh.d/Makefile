PROG =  xyc2msh
SRCS =	main.c++ \
	rotate_right.c++ \
	rotate_left.c++ \
	N_set.c++ \
	degeneracy.c++ \
	super_triangle.c++ \
	LawsonSwap.c++ \
	vanish_boundary_triangle.c++ \
	delaunay.c++ \
	count_nodes.c++ \
	generate_fN.c++ \
	fprintmesh.c++ \
	ifs2xyc.c++

OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -Wall
LDFLAGS = 

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS)
	c++ $(OBJS) -o $@ $(LDFLAGS)

.c++.o:
	c++ $(CPPFLAGS) -c $< -o $@

clean :
	rm -f $(OBJS) $(PROG) kanto.msh1 pt.msh1 test

test : $(PROG)
	./xyc2msh kanto.xyc > kanto.msh1
	./xyc2msh pt.xyc > pt.msh1
	diff kanto.msh0 kanto.msh1
	diff pt.msh0 pt.msh1
	touch test
