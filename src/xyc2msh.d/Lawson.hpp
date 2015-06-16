#ifndef LAWSON_HPP_
#define LAWSON_HPP_

#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"
#include "fplane.hpp"
#include "notLawson.hpp"

template<class Real>
int Lawson(vector<xyc>&Z, vector<nde>&N,int e, Real x, Real y)
{ 
  int n, count = 0;
  int a,b,c; Real x0,y0,x1,y1,x2,y2;
  n = N.size()-1;
  while(1){
    a = N[e].a; b = N[e].b; c = N[e].c; 
    x0 = Z[a].x; x1 = Z[b].x; x2 = Z[c].x;
    y0 = Z[a].y; y1 = Z[b].y; y2 = Z[c].y;
    if(n<count++) return notLawson(Z,N,x,y);
    
    if     (fplane(x,y,0.0, x0,y0,1.0, x1,y1,0.0, x2,y2,0.0)>0.0) e = N[e].A;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,1.0, x2,y2,0.0)>0.0) e = N[e].B;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,0.0, x2,y2,1.0)>0.0) e = N[e].C;
    else return e;
  }
}
#endif
