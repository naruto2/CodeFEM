#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "est/op.hpp"
#include "est/xmesh.hpp"


std::vector<xyc> cavity(int n)
{
  std::vector<xyc>Z;
  xyc z;

  z.x = 0.0; z.y =0.0; z.label = NULL; Z.push_back(z);

  double dx=0.5, dy=0.5;
  
  if ( n!=0 ) {
    dx = 1.0/(double)n;
    dy = 1.0/(double)n;
  }

  for ( double x = 0.0; x <1.0+dx; x+=dx) 
    for ( double y = 0.0; y<1.0+dx; y+=dy) {
      z.x = x; z.y = y; 
      z.label=NULL; 

      if ( 1.0-dy/2.0 < y && y < 1.0+dy/2.0 )
	z.label =strdup("e2");
      if ( 1.0-dx/2.0 < x && x < 1.0+dx/2.0 ) 
	z.label =strdup("e1");
      if ( y == 0.0 ) 
	z.label = strdup("e0");
      if ( x == 0.0 ) 
	z.label = strdup("e3");
      if ( x == 0.0 && y == 0.0 ) 
	z.label= strdup("v0");
      if ( x == 0.0 && 1.0-dy/2.0<y && y<1.0-dy/2.0) 
	z.label = strdup("v3");
      if ( y == 0.0 && 1.0-dx/2.0<x && x<1.0-dx/2.0)
	z.label = strdup("v1");
      if ( 1.0-dy/2.0<y && y<1.0-dy/2.0 && 1.0-dx/2.0<x && x<1.0-dx/2.0)
	z.label = strdup("v2");

      Z.push_back(z);
    }
  z.x = dx/2.0; z.y = dy/2.0; z.label = NULL; Z.push_back(z);
  z.x = 1.0-dx/2.0; z.y = 1.0-dy/2.0; z.label = NULL; Z.push_back(z);

  return Z;
}
