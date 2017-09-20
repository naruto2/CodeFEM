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
  
  M.clear();
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


void makeAx(matrix<double>&Ax,vector<double>&U,vector<xyc>&Z,vector<nde>&N){
  double del, B1, B2, B3, u[7];
  int e, m, n, i, j, I, J, a, b, c, A, B, C;
  n = N.size();
  m = dimp2(N);

  Ax.clear();
  Ax.resize(m+1);
  
  for (e=1; e<n; e++) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;

    del = delta(e,Z,N);
    i = 0;
    foreach(I, &a, &b, &c, &A, &B, &C){
      ++i; j=0;
      foreach(J, &a, &b, &c, &A, &B, &C) {
	B1 = Z[b].y - Z[c].y; B1 /= 2.0*del;
	B2 = Z[c].y - Z[a].y; B2 /= 2.0*del;
	B3 = Z[a].y - Z[b].y; B3 /= 2.0*del;
	u[1] = U[a];
	u[2] = U[b];
	u[3] = U[c];
	u[4] = U[A];
	u[5] = U[B];
	u[6] = U[C];
	Ax[I][J] += del*axij(i,++j, u, B1, B2, B3);
      }
    }
  }
}


void makeAy(matrix<double>&Ay,vector<double>&V,vector<xyc>&Z,vector<nde>&N){
  double del, C1, C2, C3, v[7];
  int e, m, n, i, j, I, J, a, b, c, A, B, C;
  n = N.size();
  m = dimp2(N);

  Ay.clear();
  Ay.resize(m+1);
  
  for (e=1; e<n; e++) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    del = delta(e,Z,N);
    i = 0;
    foreach(I, &a, &b, &c, &A, &B, &C) {
      ++i; j=0;
      foreach(J, &a, &b, &c, &A, &B, &C) {
        C1 = Z[c].x - Z[b].x; C1 /= 2.0*del;
        C2 = Z[a].x - Z[c].x; C2 /= 2.0*del;
        C3 = Z[b].x - Z[a].x; C3 /= 2.0*del;
        v[1] = V[a+m];
        v[2] = V[b+m];
        v[3] = V[c+m];
        v[4] = V[A+m];
        v[5] = V[B+m];
        v[6] = V[C+m];
        Ay[I][J] += del*ayij(i,++j, v, C1, C2, C3);
      }
    }
  }
}

int main(){
  matrix<double> M, Ax, Ay;
  vector<double> U;
  vector<xyc>Z; vector<nde>N;
  
  f2mesh(fopen("kanto.mesh","r"),Z,N);
  int i, m = dimp2(N);
  
  U.resize(2*m+Z.size());
  for (i=1; i<=2*m; i++) U[i] = 1.0;
  makeM(M,Z,N);
  makeAx(Ax,U,Z,N);
  makeAy(Ay,U,Z,N);
  plotmatrix(Ay);
  
  dij(1,1);
  hxij(1,1);
  hyij(1,1);
  return 0;
}
