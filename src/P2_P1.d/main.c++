#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "P2_P1.hpp"
#include "mij.hpp"
#include "axij.hpp"
#include "ayij.hpp"
#include "dij.hpp"
#include "hxij.hpp"
#include "hyij.hpp"
#include "est/sparse.hpp"
#include "est/xmesh.hpp"
#include "est/foreach.hpp"

using namespace sparse;

void makeM(matrix<double>&M,vector<xyc>&Z,vector<nde>&N)
{
  int e, m=dimp2(N), n=N.size(), i, j, I, J;
  int a, b, c, A, B, C;
  double s;
  
  M.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    s=delta(e,Z,N);
    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C){
      j=1;
      foreach(J,&a,&b,&c,&A,&B,&C){
	M[I][J] += mij(i,j,s);
	j++;
      }
      i++;
    }
  }
}

int main(){
  matrix<double> M;
  vector<xyc>Z; vector<nde>N;
  
  f2mesh(fopen("kanto.mesh","r"),Z,N);

  makeM(M,Z,N);
  plotmatrix(M);
  




  double u[7],v[7];
  axij(1,1,u,1.0,1.0,1.0);
  ayij(1,1,v,1.0,1.0,1.0);
  dij(1,1);
  hxij(1,1);
  hyij(1,1);
  return 0;
}
