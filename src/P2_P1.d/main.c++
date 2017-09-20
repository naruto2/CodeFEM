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


void makeD(matrix<double>&D,vector<xyc>&Z,vector<nde>&N)
{
  int e, m, n, i, j, I, J;
  int a, b, c, A, B, C;
  double del;
  
  m = dimp2(N);
  n = N.size();

  D.clear();
  D.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    del = delta(e,Z,N);
    setB1B2C1C2((Z[b].y-Z[c].y)/(2.0*del), (Z[c].y-Z[a].y)/(2.0*del),
		(Z[c].x-Z[b].x)/(2.0*del), (Z[a].x-Z[c].x)/(2.0*del));

    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C) {
      j=1;
      foreach(J,&a,&b,&c,&A,&B,&C) {
        D[I][J] += del*dij(i,j);
        j++;
      }
      i++;
    }
  }
}


void makeHx(matrix<double>&Hx,vector<xyc>&Z,vector<nde>&N)
{
  int e, m, n, i, j, I, J;
  int a, b, c, A, B, C;
  double del;
  
  m = dimp2(N);
  n = N.size();

  Hx.clear();
  Hx.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    del = delta(e,Z,N);
    setB1B2C1C2((Z[b].y-Z[c].y)/(2.0*del), (Z[c].y-Z[a].y)/(2.0*del),
		(Z[c].x-Z[b].x)/(2.0*del), (Z[a].x-Z[c].x)/(2.0*del));

    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C) {
      j=1;
      foreach(J,&a,&b,&c) {
        Hx[I][J] += del*hxij(i,j);
        j++;
      }
      i++;
    }
  }
}


void makeHy(matrix<double>&Hy,vector<xyc>&Z,vector<nde>&N)
{
  int e, m, n, i, j, I, J;
  int a, b, c, A, B, C;
  double del;
  
  m = dimp2(N);
  n = N.size();

  Hy.clear();
  Hy.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    del = delta(e,Z,N);
    setB1B2C1C2((Z[b].y-Z[c].y)/(2.0*del), (Z[c].y-Z[a].y)/(2.0*del),
		(Z[c].x-Z[b].x)/(2.0*del), (Z[a].x-Z[c].x)/(2.0*del));

    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C) {
      j=1;
      foreach(J,&a,&b,&c) {
        Hy[I][J] += del*hyij(i,j);
        j++;
      }
      i++;
    }
  }
}

double tau(void) {
  return 0.01;
}

double Re(void) {
  return 10.0;
}


void makeA(matrix<double>&A,vector<double>&U,vector<xyc>&Z,vector<nde>&N)
{
  matrix<double> M, Ax, Ay, D, Hx, Hy;
  int m;
  
  makeM(M,Z,N);
  makeAx(Ax,U,Z,N);
  makeAy(Ay,U,Z,N);
  makeD(D,Z,N);
  makeHx(Hx,Z,N);
  makeHy(Hy,Z,N);

  A.clear();
  m = dimp2(N);
  A.resize(2*m+N.size());
  
  int i, j;

  for (i=1; i<=m; i++) for (auto it : M[i]) { j = it.first;
      A[  i][  j] = M[i][j]/tau();
      A[m+i][m+j] = M[i][j]/tau();
    }

  for (i=1; i<=m; i++) for (auto it : Ax[i]) { j = it.first;
      A[  i][  j] += Ax[i][j];
      A[m+i][m+j] += Ax[i][j];
    }

  for (i=1; i<=m; i++) for (auto it : Ay[i]) { j = it.first;
      A[  i][  j] += Ay[i][j];
      A[m+i][m+j] += Ay[i][j];
    }

  for (i=1; i<=m; i++) for (auto it : D[i]) { j = it.first;
      A[  i][  j] += D[i][j]/Re();
      A[m+i][m+j] += D[i][j]/Re();
    }

  for (i=1; i<=m; i++) for (auto it : Hx[i]) { j = it.first;
      A[    i][2*m+j] = -Hx[i][j];
      A[2*m+j][    i] = -Hx[i][j];
    }

  for (i=1; i<=m; i++) for (auto it : Hy[i]) { j = it.first;
      A[  m+i][2*m+j] = -Hy[i][j];
      A[2*m+j][  m+i] = -Hy[i][j];
    }
}


int main(){
  matrix<double> A;
  vector<double> U;
  vector<xyc>Z; vector<nde>N;
  
  f2mesh(fopen("cavity.mesh","r"),Z,N);
  int i, m = dimp2(N);
  
  U.resize(2*m+Z.size());
  for (i=1; i<=2*m; i++) U[i] = 1.0;

  makeA(A,U,Z,N);

  system("date");
  plotmatrix(A);
  return 0;
}
