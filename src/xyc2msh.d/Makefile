PROG =  xyc2msh
SRCS =	main.c++ \
	distance2.c++ \
	sarrus.c++ \
	fplane.c++ \
	cramer3.c++ \
	rotate_right.c++ \
	rotate_left.c++ \
	incircle.c++ \
	N_set.c++ \
	notLawson.c++ \
	Lawson.c++ \
	degeneracy.c++ \
	min2.c++ \
	max2.c++ \
	super_triangle.c++ \
	LawsonSwap.c++ \
	vanish_boundary_triangle.c++ \
	delaunay1.c++ \
	count_nodes.c++ \
	generate_fN.c++ \
	arytovec.c++ \
	fprintmesh.c++ \
	fp2xyc.c++

OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS = -Wall -g -fpermissive
LDFLAGS = -lestiva -lm

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
