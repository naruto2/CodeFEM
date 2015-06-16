#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"

void vanish_boundary_triangle(vector<xyc>&Z, vector<nde>&N)
{

  int A,B,C,i;

  for (i=1; i<(int)N.size(); i++) {
    if(Z[N[i].a].label==NULL||Z[N[i].b].label==NULL||Z[N[i].c].label==NULL)
      continue;

    
    A = N[i].A;
    B = N[i].b;
    C = N[i].c;
    if(N[A].A == i) N[A].A=0;
    if(N[A].B == i) N[A].B=0;
    if(N[A].C == i) N[A].C=0;
    if(N[B].A == i) N[B].A=0;
    if(N[B].B == i) N[B].B=0;
    if(N[B].C == i) N[B].C=0;
    if(N[C].A == i) N[C].A=0;
    if(N[C].B == i) N[C].B=0;
    if(N[C].C == i) N[C].C=0;

    N[i].a=0;
    N[i].b=0;
    N[i].c=0;
    N[i].A=0;
    N[i].B=0;
    N[i].C=0;
  }
}
