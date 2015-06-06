#include <vector>
using namespace std;
#include "xyc_nde.h"

#define T double
extern T distance2(T,T,T,T);
extern void cramer3(T*,T*,T*,T,T,T,T,T,T,T,T,T,T,T,T);

int incircle(vector<xyc>&Z, vector<nde> &N, int p,int e2)
{
  T a, b, c, x,y,x0,y0,x1,y1,x2,y2;

  x=Z[p].x; y=Z[p].y;
  x0=Z[N[e2].a].x;  y0=Z[N[e2].a].y;
  x1=Z[N[e2].b].x;  y1=Z[N[e2].b].y;
  x2=Z[N[e2].c].x;  y2=Z[N[e2].c].y;

  cramer3(&a,&b,&c, 
	  x0,y0,1.0,
	  x1,y1,1.0,
	  x2,y2,1.0, 
	  x0*x0 + y0*y0,
	  x1*x1 + y1*y1, 
	  x2*x2 + y2*y2);

  return distance2(x0,y0,a/2.0,b/2.0)-distance2(x,y,a/2.0,b/2.0)>0.0? 1:0;

}
