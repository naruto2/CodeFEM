#ifndef _EST_XYC2MSH_HPP_
#define _EST_XYC2MSH_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

using namespace std;

typedef struct{ int a, b, c, A, B, C; } nde;
typedef struct{ double x, y; char *label; } xyc;

extern void in2xyc( istream &, vector<xyc>& );
extern void delaunay( vector<xyc>&, vector<nde>& );
extern void outmesh( ostream &, vector<xyc>&, vector<nde>& );
extern void plotmesh( vector<xyc>&, vector<nde>& );
extern void sortmesh( vector<xyc>&, vector<nde>& );
extern void p2( vector<xyc>&, vector<nde>& );
extern void f2mesh(FILE*, vector<xyc>&, vector<nde>&);

#endif
