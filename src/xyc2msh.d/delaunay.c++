#include <cstdio>
#include <vector>
using namespace std;
#include "xyc_nde.h"

extern void super_triangle(vector<xyc>&);
extern void LawsonSwap(vector<xyc>&, vector<nde>&);
extern void vanish_boundary_triangle(vector<xyc>&, vector<nde>&);
extern int count_nodes(vector<nde>&);
extern int generate_fN(vector<nde>&, vector<int>&);

void delaunay(vector<xyc>&Zo, vector<nde>&No)
{ int i, n, z; 
  vector<xyc> Z;
  vector<nde> N, N0;
  vector<int> fN, fNinv, fZ,fZinv;
  
  z = Zo.size()-1;
  Z.resize(z+1);
  for ( i=0; i<z; i++) Z[i] = Zo[i];
  
  super_triangle(Z);
  N.resize(Z.size()*3);
  LawsonSwap(Z,N);
  vanish_boundary_triangle(Z, N);
  
  fN.resize(count_nodes(N)+1);
  n = generate_fN(N,fN);
  
  N0.resize(fN.size());
  for (i=1;i<(int)N0.size(); i++) N0[i]=N[fN[i]];
  
  fNinv.resize(n);
  fZ.resize(z+1);

  for (i=1; i<(int)fN.size(); i++) fNinv[fN[i]]=i;

  for (i=1; i<=(int)N0.size(); i++) {
    N0[i].A = fNinv[N0[i].A];
    N0[i].B = fNinv[N0[i].B];
    N0[i].C = fNinv[N0[i].C];
  }  

  for(z=0, i=1; i<(int)N0.size(); i++) {
    if (fZ[N0[i].a] == 0) fZ[N0[i].a] = ++z;
    if (fZ[N0[i].b] == 0) fZ[N0[i].b] = ++z;
    if (fZ[N0[i].c] == 0) fZ[N0[i].c] = ++z;
  }
    
  for(i=1; i<(int)N0.size(); i++) {
    N0[i].a = fZ[N0[i].a];
    N0[i].b = fZ[N0[i].b];
    N0[i].c = fZ[N0[i].c];
  }  

  fZinv.resize(z+1);
  for (i=1; i<(int)fZ.size(); i++) fZinv[fZ[i]] = i;

  Zo.resize(z+2);
  for (i=1; i<=z; i++) Zo[i]  = Z[fZinv[i]];  

  No.resize(N0.size()+1);
  for (i=0;i<(int)N0.size();i++) No[i] = N0[i];
  
}

