#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
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

void printmx(sparse::matrix<double>&A){
  for (int i=1; i<A.size(); i++) for ( auto it: A[i]) { int j = it.first;
      printf("%d %d %e\n",i,j,A[i][j]);
    }
  exit(0);
}

void makeM(sparse::matrix<double>&M,vector<xyc>&Z,vector<nde>&N)
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


void makeAx(sparse::matrix<double>&Ax,vector<double>&U,vector<xyc>&Z,vector<nde>&N){
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


void makeAy(sparse::matrix<double>&Ay,vector<double>&V,vector<xyc>&Z,vector<nde>&N){
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


void makeD(sparse::matrix<double>&D,vector<xyc>&Z,vector<nde>&N)
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


void makeHx(sparse::matrix<double>&Hx,vector<xyc>&Z,vector<nde>&N)
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


void makeHy(sparse::matrix<double>&Hy,vector<xyc>&Z,vector<nde>&N)
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
  return 0.07 ;
}

double Re(void) {
  return 1000.0;
}

static char* border(char *s, char *t)
{
  return (strcmp(s,t)<0?s:t);
}

void makeMid(vector<xyc>&Mid,vector<xyc>&Z,vector<nde>&N) {
  int e, n=N.size();

  Mid.clear();
  Mid.resize(dimp2(N)+1);
  
  for(e=1; e<n; e++) {
    Mid[N[e].A].x = (Z[N[e].b].x + Z[N[e].c].x)/2.0;
    Mid[N[e].A].y = (Z[N[e].b].y + Z[N[e].c].y)/2.0;
    Mid[N[e].A].label = border(Z[N[e].b].label,Z[N[e].c].label);

    Mid[N[e].B].x = (Z[N[e].c].x + Z[N[e].a].x)/2.0;
    Mid[N[e].B].y = (Z[N[e].c].y + Z[N[e].a].y)/2.0;
    Mid[N[e].B].label = border(Z[N[e].c].label,Z[N[e].a].label);

    Mid[N[e].C].x = (Z[N[e].a].x + Z[N[e].b].x)/2.0;
    Mid[N[e].C].y = (Z[N[e].a].y + Z[N[e].b].y)/2.0;
    Mid[N[e].C].label = border(Z[N[e].a].label,Z[N[e].b].label);
  }

  for(int i=1;i<Z.size();i++){
    Mid[i].x = Z[i].x;
    Mid[i].y = Z[i].y;
    Mid[i].label = Z[i].label;
  }
}  

void makeA(sparse::matrix<double>&A,vector<double>&U,vector<double>&b,vector<xyc>&Z,vector<nde>&N,vector<xyc>&Mid)
{
  sparse::matrix<double> M, Ax, Ay, D, Hx, Hy;
  int m;
  
  makeM(M,Z,N);
  makeAx(Ax,U,Z,N);
  makeAy(Ay,U,Z,N);
  makeD(D,Z,N);
  makeHx(Hx,Z,N);
  makeHy(Hy,Z,N);
  
  A.clear();
  m = dimp2(N);
  A.resize(2*m+Z.size());
  
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
  
  b.clear();

  for (i=1; i<=m; i++) {
    b[  i] = 0.0;
    b[m+i] = 0.0;
    for (j=1; j<=m; j++) { 
      b[  i] += M[i][j]*U[  j]/tau();
      b[m+i] += M[i][j]*U[m+j]/tau();
    }
  }
  
  m = dimp2(N);
  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v0")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v1")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v2")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v3")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e0")) {
      A[i+m].clear();
      A[i+m][i+m] = 1.0;
      b[i+m] = 0.0;
    }


  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e1")) {
      A[i].clear();
      A[i][i] = 1.0;
      b[i] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e2")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 1.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e3")) {
      A[i].clear();
      A[i][i] = 1.0;
      b[i] = 0.0;
    }

  A[2*m+1].clear();
  A[2*m+1][2*m+1] = 1.0;
  b[2*m+1] = 0.0;
  
}

void plotuv(FILE *pp,vector<double>&U,vector<xyc>&Z,vector<nde>&N,vector<xyc>&Mid){

  double scale=0.4;
  long arrow =1 ;
  int i, m=Mid.size()-1;
  
  for (i=1; i<=m; i++)
    fprintf(pp, "set arrow %ld from %f,%f to %f,%f\n",
	    arrow++,Mid[i].x,Mid[i].y,Mid[i].x+U[i]*scale,Mid[i].y+U[i+m]*scale);
  fprintf(pp,"set xrange [0:1]\n");
  fprintf(pp,"set yrange [0:1]\n");
  fprintf(pp,"plot '-' w l\n");
  for(int e=1;e<N.size();e++){
    fprintf(pp,"%f %f\n",Z[N[e].a].x,Z[N[e].a].y);
    fprintf(pp,"%f %f\n",Z[N[e].b].x,Z[N[e].b].y);
    fprintf(pp,"%f %f\n",Z[N[e].c].x,Z[N[e].c].y);
    fprintf(pp,"%f %f\n\n",Z[N[e].a].x,Z[N[e].a].y);
  }
  fprintf(pp,"e\n");
  fflush(pp);
  sleep(1);
}

#include "est/solver.hpp"
#include "est/matrix.hpp"


void sparse__solve(sparse::matrix<double>&A,vector<double>&U,vector<double>&b)
{
  int i, j;
  matrix AA(A.size()-1);
  vector<double> x(A.size()-1), bb(A.size()-1);
  Preconditioner M;
      
  for ( i=1; i<A.size(); i++) for (auto it: A[i]) { j = it.first;
	AA[i-1][j-1] = A[i][j];
      }
    for ( i=1; i<A.size(); i++)
      bb[i-1] = b[i];

    x = bicgstab(M,AA,bb);
    vector<double> gpubicgstab(matrix&, vector<double>&);
    //x = gpubicgstab(AA,bb);

    for(i=0; i<=x.size(); i++){
      U[i+1] = x[i];
    }

}


int main(){
  FILE *pp = popen("gnuplot","w");
  
  vector<xyc>Z; vector<nde>N;
  f2mesh(fopen("cavity16.mesh","r"),Z,N);

  vector<xyc> Mid;
  makeMid(Mid,Z,N);

  int i, j, m = dimp2(N);

  sparse::matrix<double> A;
  A.clear();
  A.resize(2*m+Z.size());
  vector<double> U;
  U.clear();
  U.resize(2*m+Z.size());
  vector<double> b;
  b.clear();
  b.resize(2*m+Z.size());


  for(int k=0;k<1000;k++){
    fprintf(stderr,"k=%d\n",k);
    makeA(A,U,b,Z,N,Mid);
    sparse__solve(A,U,b);
    plotuv(pp,U,Z,N,Mid);
  }
  sleep(30);
}