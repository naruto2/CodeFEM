#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"
#include "ary.h"


extern void arytovec_xyc(xyc*, vector<xyc>&);
extern void arytovec_nde(nde*, vector<nde>&);
extern void vectoary_xyc(vector<xyc> &, xyc *);
extern void vectoary_nde(vector<nde> &, nde *);

extern void super_triangle(vector<xyc>&);
extern void LawsonSwap(vector<xyc>&, vector<nde>&);
extern void vanish_boundary_triangle(vector<xyc>&, vector<nde>&);
extern int count_nodes(vector<nde>&);
extern int generate_fN(vector<nde>&,int*);

void delaunay1(xyc **Zo, nde **No)
{ int i, n, z; 
  static int *fN,*fNinv,*fZ,*fZinv;
  static xyc *Z0, *Z;
  static nde *N0, *N;
  vector<xyc> Zv;
  vector<nde> Nv;
  
  Z = *Zo;
  
  arytovec_xyc(Z,Zv);
  super_triangle(Zv);

  ary1(N,(Zv.size()-1)*3);
  arytovec_nde(N,Nv);

  LawsonSwap(Zv,Nv);
  vanish_boundary_triangle(Zv, Nv);

  ary1(fN,count_nodes(Nv)+1);



  n = generate_fN(Nv,fN);

  vectoary_xyc(Zv,Z);
  vectoary_nde(Nv,N);

  
  ary1(N0, dim1(fN)+1); for (i=1;i<=dim1(N0); i++) N0[i]=N[fN[i]];
  
  ary1(fNinv,n+1); 
  for (i=1; i<=dim1(fN); i++) fNinv[fN[i]]=i;

  for (i=1; i<=dim1(N0); i++) {
    N0[i].A = fNinv[N0[i].A];
    N0[i].B = fNinv[N0[i].B];
    N0[i].C = fNinv[N0[i].C];
  }  

  ary1(fZ,dim1(Z)+1); 

  for(z=0, i=1; i<=dim1(N0); i++) {
    if (fZ[N0[i].a] == 0) fZ[N0[i].a] = ++z;
    if (fZ[N0[i].b] == 0) fZ[N0[i].b] = ++z;
    if (fZ[N0[i].c] == 0) fZ[N0[i].c] = ++z;
  }
    
  for(i=1; i<=dim1(N0); i++) {
    N0[i].a = fZ[N0[i].a];
    N0[i].b = fZ[N0[i].b];
    N0[i].c = fZ[N0[i].c];
  }  

  ary1(fZinv,z+1); for (i=1; i<=dim1(fZ); i++) fZinv[fZ[i]] = i;
  ary1(Z0,z+1);    for (i=1; i<=dim1(Z0); i++) Z0[i] = Z[fZinv[i]];
    
  *Zo = Z0;
  *No = N0;
}
