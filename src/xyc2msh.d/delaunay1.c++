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
extern int generate_fN(vector<nde>&, vector<int>&);

void delaunay1(vector<xyc>&Zov, vector<nde>&Nov)
{ int i, n, z; 
  vector<xyc> Zv;
  vector<nde> Nv, N0v;
  vector<int> fNv, fNinvv, fZv,fZinvv;
  
  z = Zov.size()-1;
  Zv.resize(z+1);
  for ( i=0; i<z; i++) Zv[i] = Zov[i];
  
  super_triangle(Zv);
  Nv.resize(Zv.size()*3);
  LawsonSwap(Zv,Nv);
  vanish_boundary_triangle(Zv, Nv);
  
  fNv.resize(count_nodes(Nv)+1);
  n = generate_fN(Nv,fNv);
  
  N0v.resize(fNv.size());
  for (i=1;i<(int)N0v.size(); i++) N0v[i]=Nv[fNv[i]];
  
  fNinvv.resize(n);
  fZv.resize(z+1);

  for (i=1; i<(int)fNv.size(); i++) fNinvv[fNv[i]]=i;

  for (i=1; i<=(int)N0v.size(); i++) {
    N0v[i].A = fNinvv[N0v[i].A];
    N0v[i].B = fNinvv[N0v[i].B];
    N0v[i].C = fNinvv[N0v[i].C];
  }  

  for(z=0, i=1; i<(int)N0v.size(); i++) {
    if (fZv[N0v[i].a] == 0) fZv[N0v[i].a] = ++z;
    if (fZv[N0v[i].b] == 0) fZv[N0v[i].b] = ++z;
    if (fZv[N0v[i].c] == 0) fZv[N0v[i].c] = ++z;
  }
    
  for(i=1; i<(int)N0v.size(); i++) {
    N0v[i].a = fZv[N0v[i].a];
    N0v[i].b = fZv[N0v[i].b];
    N0v[i].c = fZv[N0v[i].c];
  }  

  fZinvv.resize(z+1);

  for (i=1; i<(int)fZv.size(); i++) fZinvv[fZv[i]] = i;


  Zov.resize(z+2);
  for (i=1; i<=z; i++) Zov[i]  = Zv[fZinvv[i]];  

  Nov.resize(N0v.size()+1);
  for (i=0;i<(int)N0v.size();i++) Nov[i] = N0v[i];
  
}
