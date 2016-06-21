#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include "est/op.hpp"
#include "est/xmesh.hpp"
#include "est/matrix.hpp"
#include "est/solver.hpp"

using namespace std;

extern vector<xyc> ncpolynomial1(vector<xyc> Z, vector<nde> &N );
extern void estiva_forgammap1(long *x);
extern int estiva_forgammap1_loop(long *x, const char *name, vector<xyc> &Z);
extern void printmatrix(matrix &A, const char *name);  
extern void printvector(vector<double> &b, const char *name);
extern void plotncpolynomial1(FILE *fp, vector<xyc> Mid, vector<double> u, vector<double> v);

#define length(a,b) \
  ((Z[b].x-Z[a].x)*(Z[b].x-Z[a].x)+(Z[b].y-Z[a].y)*(Z[b].y-Z[a].y))

#define inner(a,b,c) \
  ((Z[c].x-Z[a].x)*(Z[b].x-Z[c].x)+(Z[c].y-Z[a].y)*(Z[b].y-Z[c].y))


static vector<double> S_(vector<xyc> Z, vector<nde> N)
{

  long i, e;
  double xi,xj,xk,yi,yj,yk;

  e = N.size()-1;
  vector<double> S;
  S.resize(e+1);
  
  for(i=1;i<=e;i++){
    xi = Z[N[i].a].x;
    xj = Z[N[i].b].x;
    xk = Z[N[i].c].x;
    yi = Z[N[i].a].y;
    yj = Z[N[i].b].y;
    yk = Z[N[i].c].y;
    S[i] = xi*yj+xj*yk+xk*yi-yi*xj-yj*xk-yk*xi;
    S[i] /= 2.0;
  }
  return S;
}


static matrix M__(vector<xyc> Mid, vector<nde> N, vector<double> S)
{
  matrix M;
  long  e, a, b, c, A, B, C, m, n;

  m = Mid.size()-1;
  n = N.size()-1;
  M.resize(m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    M[A][A] += S[e]/3.0;
    M[B][B] += S[e]/3.0;
    M[C][C] += S[e]/3.0;
  }
  return M;
}

static matrix K__(vector<xyc> Mid, vector<xyc> Z, vector<nde> N, vector<double> S)
{
  matrix K;
  long e, a, b, c, A, B, C, m, n;
  double s;

  m = Mid.size()-1;
  n = N.size()-1;
  K.resize(m+1);
  
  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    s = S[e];

    K[A][A] += length(b,c)/s; K[A][B] +=inner(a,b,c)/s; K[A][C] +=inner(a,c,b)/s;
    K[B][A] +=inner(b,a,c)/s; K[B][B] += length(c,a)/s; K[B][C] +=inner(b,c,a)/s;
    K[C][A] +=inner(c,a,b)/s; K[C][B] +=inner(c,b,a)/s; K[C][C] += length(a,b)/s;

  }
  return K;
}

static matrix Hx__(vector<xyc> Mid, vector<xyc> Z, vector<nde> N)
{
  matrix Hx;
  long e, a, b, c, A, B, C, m, n;

  m = Mid.size()-1;
  n = N.size()-1;
  Hx.resize(m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hx[A][e] += -(Z[c].y-Z[b].y);
    Hx[B][e] += -(Z[a].y-Z[c].y);
    Hx[C][e] += -(Z[b].y-Z[a].y);
  }
  return Hx;
}


static matrix Hy__(vector<xyc> Mid, vector<xyc> Z, vector<nde> N)
{
  matrix Hy;
  long e, a, b, c, A, B, C, m, n;

  m = Mid.size()-1;
  n = N.size()-1;
  Hy.resize(m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hy[A][e] += (Z[c].x-Z[b].x);
    Hy[B][e] += (Z[a].x-Z[c].x);
    Hy[C][e] += (Z[b].x-Z[a].x);

  }
  return Hy;
}


static void A__(matrix &A, vector<xyc> &Mid, vector<nde> &N, matrix &M, double tau, matrix &K, matrix &Hx, matrix &Hy)
{
  long   i, j, NUM, m, n;

  m   = Mid.size()-1;
  n   = N.size()-1;
  NUM = m*2+n;
  A.resize(NUM+1);

  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      A[  i][  j] = M[i][j] + tau*K[i][j];
      A[m+i][m+j] = M[i][j] + tau*K[i][j];
    }
  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= n; j++ ) {
      A[    i][2*m+j] = -tau*Hx[i][j];
      A[2*m+j][    i] = -tau*Hx[i][j];
      A[  m+i][2*m+j] = -tau*Hy[i][j];
      A[2*m+j][  m+i] = -tau*Hy[i][j];
    }
}

void Rhs(vector<double> &b, vector<xyc> &Mid, vector<nde> &N,matrix &M,double t,vector<double> &Fx,vector<double> &Fy,
	 vector<double> &Ux, vector<double> &Uy, vector<double> &x)
{
  long   i, j, NUM, m, n;

  m = Mid.size()-1;
  n = N.size()-1;

  NUM = m*2+n;

  for( i = 1; i <= NUM; i++ ) b[i] = 0.0;

  for( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      b[  i] += M[i][j]*x[  j];
      b[m+i] += M[i][j]*x[m+j];
    }

}

#define forgammap1(i,NAME,Z)                                            \
  for ( estiva_forgammap1(&(i));      estiva_forgammap1_loop(&(i),NAME,Z);)


static void zerofillrow(matrix &A, unsigned long i)
{
  unsigned long j;

  for ( j = 0; j < A.size(); j++ ) {
    A[i][j] = 0.0;
  }
}

static void boundary_condition(vector<nde> &N, vector<xyc> &Mid, matrix &A, vector<double> &b)
{
  long NUM, i, j, m, n;
  m = Mid.size()-1;
  n = N.size()-1;
  NUM = 2*m+n;

  printf ("\n---------- e0 ------------\n");
  forgammap1(i,"e0",Mid){
    zerofillrow(A,i+m);
    A[i+m][i+m] = 1.0;
    b[i+m] = 0.0;
    printf("%ld ",i);
  }

  printf ("\n---------- e1 ------------\n");
  forgammap1(i,"e1",Mid){
    zerofillrow(A,i);
    A[i][i] = 1.0;
    b[i] = 0.0;
    printf("%ld ",i);
  }

  printf ("\n---------- e2 ------------\n");

  forgammap1(i,"e2",Mid){
    zerofillrow(A,i);
    A[i][i] = 1.0;
    b[i]   = 0.1;
    zerofillrow(A,i+m);
    A[i+m][i+m] = 1.0;
    b[i+m]     = 0.0;
    printf("%ld ",i);
  }

  printf ("\n---------- e3 ------------\n");

  forgammap1(i,"e3",Mid){
    zerofillrow(A,i);
    A[i][i] = 1.0;
    b[i] = 0.0;
    printf("%ld ",i);
  }

  printf ("\n--------------------------\n");

  zerofillrow(A,NUM-1);
  A[NUM-1][NUM-1] = 1.0;
  b[NUM-1] = 0.0001;
}



int main(int argc, char ** argv)
{
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  ifstream ifs(argv[1]);

  if (!ifs) {
    cerr << "Error: 入力ストリームを開けませんでした(xmesh)" << endl;
    return 0;
  }

  in2xyc(ifs,Z);
  ifs.close();

  delaunay(Z, N);
  printf("N.size()-1 = %ld\n",N.size()-1);
  
  vector<xyc> Mid = ncpolynomial1(Z,N);

  for ( unsigned long i=1; i< Mid.size(); i++ ) {
    printf("%d %f %f %s\n",i, Mid[i].x, Mid[i].y, Mid[i].label );
  }
  printf("-------------------------------------\n");

  for ( unsigned long e=1; e<N.size(); e++ ) {
    printf("%d %d %d %d %d %d %d\n",e,N[e].a,N[e].b,N[e].c,N[e].A,N[e].B,N[e].C);
  }    
  printf("-------------- N -----------------------\n");
  
  vector<double> S = S_(Z,N);
  
  for ( unsigned long e=1; e<N.size(); e++ ) {
    printf("%ld %f\n",e,S[e]);
  }

  unsigned long m = Mid.size()-1;
  unsigned long n = S.size()-1;
  unsigned long NUM = 2*m+n;

  matrix M = M__(Mid, N, S);
  matrix K = K__(Mid, Z, N, S);
  matrix Hx= Hx__(Mid, Z, N);
  matrix Hy= Hy__(Mid, Z, N);
  double t = 0.001;
  matrix A;
  vector<double> Fx(m+1), Fy(m+1), Ux(m+1), Uy(m+1), x(NUM+1), b(NUM+1);
  
  unsigned long i, k;
  for ( k = 1; k<=1000; k++ ) {
    A__(A, Mid, N, M,t,K,Hx,Hy);
    Rhs(b, Mid, N, M, t, Fx, Fy, Ux, Uy, x);
    boundary_condition(N,Mid,A,b);
    A[0][0] = 1.0;
    x = solve(A,b);
    for(i=1;i<=m;i++){ Ux[i] = x[i]; Uy[i] = x[i+m];}
    printf("k = %ld\n",k);
    if ( k == 1000 ) {
      plotncpolynomial1(fopen("foo","w"), Mid, Ux, Uy);
      exit(0);
    }
  }
  return 0;
}
