#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"

#define T double
extern T fplane(T,T,T,T,T,T,T,T,T,T,T,T);
extern int notLawson(vector<xyc>&,vector<nde>&,T,T);
  
int Lawson(vector<xyc>&Z, vector<nde>&N,int e, T x, T y)
{ 
  int n, count = 0;
  int a,b,c; T x0,y0,x1,y1,x2,y2;
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
