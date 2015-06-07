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
  static int *fN,*fNinv,*fZ;
  static xyc *Z0, *Z;
  static nde *N0;
  vector<xyc> Zv;
  vector<nde> Nv;
  
  Z = *Zo;
  
  arytovec_xyc(Z,Zv);
  super_triangle(Zv);

  Nv.resize(Zv.size()*3);

  LawsonSwap(Zv,Nv);
  vanish_boundary_triangle(Zv, Nv);
  
  ary1(fN,count_nodes(Nv)+1);
  vector<int> fNv, fNinvv, fZv;
  vector<nde> N0v;


  n = generate_fN(Nv,fN);
  fNv.resize(dim1(fN)+1);
  for ( i = 0; i<=dim1(fN); i++ ) fNv[i] = fN[i];

  ary1(N0, fNv.size());
  N0v.resize(fNv.size()+1);
  for (i=1;i<=dim1(N0); i++) N0[i]=Nv[fNv[i]];
  arytovec_nde(N0,N0v);
  for (i=1;i<=dim1(N0); i++) N0v[i]=Nv[fNv[i]];

  
  ary1(fNinv,n+1); 
  ary1(fZ,dim1(Z)+1);
  fNinvv.resize(n);
  fZv.resize(dim1(Z));

  for (i=0; i<(int)fNinvv.size(); i++) fNinvv[i] = fNinv[i];
  for (i=0; i<(int)fZv.size(); i++)    fZv[i] = fZ[i];
  for (i=1; i<(int)fNv.size(); i++) fNinv[fNv[i]]=i;
  for (i=1; i<(int)fNv.size(); i++) fNinvv[fNv[i]]=i;

  for (i=1; i<=(int)N0v.size(); i++) {
    N0v[i].A = fNinvv[N0v[i].A];
    N0v[i].B = fNinvv[N0v[i].B];
    N0v[i].C = fNinvv[N0v[i].C];
  }  

  for(z=0, i=1; i<(int)N0v.size(); i++) {
    if (fZ[N0v[i].a] == 0) fZ[N0v[i].a] = ++z;
    if (fZ[N0v[i].b] == 0) fZ[N0v[i].b] = ++z;
    if (fZ[N0v[i].c] == 0) fZ[N0v[i].c] = ++z;
  }
    
  for(i=1; i<(int)N0v.size(); i++) {
    N0v[i].a = fZ[N0v[i].a];
    N0v[i].b = fZ[N0v[i].b];
    N0v[i].c = fZ[N0v[i].c];
  }  

  vector<int> fZinvv;
  fZinvv.resize(z+1);

  for (i=1; i<=dim1(fZ); i++) fZinvv[fZ[i]] = i;

  ary1(Z0,z+1);
  for (i=1; i<=dim1(Z0); i++) Z0[i] = Zv[fZinvv[i]];  
  *Zo = Z0;

  for (i=0;i<(int)N0v.size();i++) N0[i] = N0v[i];
  *No = N0;
}
