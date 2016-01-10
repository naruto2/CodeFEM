#include <vector>
using namespace std;
#include "est/xmesh.hpp"
#include "distance2.hpp"

extern void rotate_right(vector<nde>&,int);
extern void rotate_left( vector<nde>&,int);

int degeneracy(vector<xyc>&Z, vector<nde>&N,int e1,int e2)
{ int a, d, W, Y; double x1,y1,x2,y2,x3,y3,x4,y4;

  if(N[e1].B == e2) rotate_left(N,e1); if(N[e1].C == e2) rotate_right(N,e1);
  if(N[e2].B == e1) rotate_left(N,e2); if(N[e2].C == e1) rotate_right(N,e2);

  x1= Z[N[e1].a].x;  y1= Z[N[e1].a].y;
  x2= Z[N[e2].a].x;  y2= Z[N[e2].a].y;
  x3= Z[N[e2].b].x;  y3= Z[N[e2].b].y;
  x4= Z[N[e2].c].x;  y4= Z[N[e2].c].y;
 
  if(distance2(x1,y1,x2,y2)==distance2(x3,y3,x4,y4)) return 1;
  
  a=N[e1].a; Y=N[e1].C;
  d=N[e2].a; W=N[e2].C;
  N[e1].b=d; N[e1].A=W; N[e1].C=e2;
  N[e2].b=a; N[e2].A=Y; N[e2].C=e1;

  if(N[W].B == e2) rotate_left(N,W); if(N[W].C == e2) rotate_right(N,W);
  if(N[Y].B == e1) rotate_left(N,Y); if(N[Y].C == e1) rotate_right(N,Y);
  N[W].A = e1; N[Y].A = e2;

  return 0;
}
