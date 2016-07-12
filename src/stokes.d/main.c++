#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <unistd.h>
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
extern void plotncpolynomial1(vector<xyc> Mid, vector<double> x);
extern void setanimefilename(const char *fname);
extern void squaremesh(int n, vector<xyc> &Z);
extern void start_baton(void);
extern int end_baton(void);
#define batonth() for(start_baton();end_baton();)

#define length(a,b) \
  ((Z[b].x-Z[a].x)*(Z[b].x-Z[a].x)+(Z[b].y-Z[a].y)*(Z[b].y-Z[a].y))

#define inner(a,b,c) \
  ((Z[c].x-Z[a].x)*(Z[b].x-Z[c].x)+(Z[c].y-Z[a].y)*(Z[b].y-Z[c].y))


static vector<double> S_(vector<xyc> &Z, vector<nde> &N)
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


static matrix M__(vector<xyc> &Mid, vector<nde> &N, vector<double> &S)
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

static matrix K__(vector<xyc> &Mid, vector<xyc> &Z, vector<nde> &N, vector<double> &S)
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

static matrix Hx__(vector<xyc> &Mid, vector<xyc> &Z, vector<nde> &N)
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


static matrix Hy__(vector<xyc> &Mid, vector<xyc> &Z, vector<nde> &N)
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
  A.clear();
  A.resize(NUM+1);

  for ( i = 1; i <= m; i++ ) for ( auto it : M[i] ) { int j = it.first;
      A[  i][  j] = M[i][j];
      A[m+i][m+j] = M[i][j];
    }
  for ( i = 1; i <= m; i++ ) for ( auto it : K[i] ) { int j = it.first;
      A[  i][  j] += tau*K[i][j];
      A[m+i][m+j] += tau*K[i][j];
    }
  for ( i = 1; i <= m; i++ ) for ( auto it : Hx[i] ) { int j = it.first;
      A[    i][2*m+j] = -tau*Hx[i][j];
      A[2*m+j][    i] = -tau*Hx[i][j];
    }
  for ( i = 1; i <= m; i++ ) for ( auto it : Hy[i] ) { int j = it.first;
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

  for( i = 1; i <= m; i++ )  for ( auto it : M[i] ) { int j = it.first;
      b[  i] += M[i][j]*x[  j];
      b[m+i] += M[i][j]*x[m+j];
    }
}

#define forgammap1(i,NAME,Z)                                            \
  for ( estiva_forgammap1(&(i)); estiva_forgammap1_loop(&(i),NAME,Z);)


static void zerofillrow(matrix &A, unsigned long i)
{
  A[i].clear();
}

static void boundary_condition(vector<nde> &N, vector<xyc> &Mid, matrix &A, vector<double> &b)
{
  long NUM, i, j, m, n;
  m = Mid.size()-1;
  n = N.size()-1;
  NUM = 2*m+n;

  forgammap1(i,"e0",Mid){
    zerofillrow(A,i+m);
    A[i+m][i+m] = 1.0;
    b[i+m] = 0.0;
  }
  forgammap1(i,"e1",Mid){
    zerofillrow(A,i);
    A[i][i] = 1.0;
    b[i] = 0.0;
  }
  forgammap1(i,"e2",Mid){
    zerofillrow(A,i);
    A[i][i] = 1.0;
    b[i]   = 0.5;
    zerofillrow(A,i+m);
    A[i+m][i+m] = 1.0;
    b[i+m]     = 0.0;
  }
  forgammap1(i,"e3",Mid){
    zerofillrow(A,i);
    A[i][i] = 1.0;
    b[i] = 0.0;
  }
  zerofillrow(A,NUM-1);
  A[NUM-1][NUM-1] = 1.0;
  b[NUM-1] = 0.0001;

  A[0][0] = 1.0;
}


#include <termios.h>
#include <unistd.h>
#include <fcntl.h>

int kbhit(void)
{
  struct termios oldt, newt;
  int ch;
  int oldf;

  tcgetattr(STDIN_FILENO, &oldt);
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

  ch = getchar();

  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  fcntl(STDIN_FILENO, F_SETFL, oldf);

  if (ch != EOF) {
    ungetc(ch, stdin);
    return 1;
  }

  return 0;
}

void mkcont(int argc, char**argv);
void contop(int &argc0, char** &argv0);

int main(int argc, char ** argv)
{
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  {
    unsigned long n = 4;
    if ( defop("-n") ) n = atoi(getop("-n").c_str());
    squaremesh(n,Z);  
  }
  delaunay(Z, N);
  sortmesh(Z,N);
  
  vector<xyc>   Mid = ncpolynomial1(Z,N);
  vector<double>  S = S_(Z,N);
  unsigned long   m = Mid.size()-1;
  unsigned long   n = S.size()-1;
  unsigned long NUM = 2*m+n;

  matrix M = M__(Mid, N, S);
  matrix K = K__(Mid, Z, N, S);
  matrix Hx= Hx__(Mid, Z, N);
  matrix Hy= Hy__(Mid, Z, N);
  double t = 0.001;
  if ( defop("-t") ) t = atof(getop("-t").c_str());
  vector<double> Fx(m+1), Fy(m+1), Ux(m+1), Uy(m+1), x(NUM+1), b(NUM+1);
  matrix A;  
  unsigned long i, k;
  setbuf(stdout,NULL);
  for ( k = 1; k<=60; k++ ) {
    A__(A, Mid, N, M,t,K,Hx,Hy);
    Rhs(b, Mid, N, M, t, Fx, Fy, Ux, Uy, x);
    boundary_condition(N,Mid,A,b);

    printf("\nk = %ld  ",k);
    /* batonth() */{

      auto f = async(launch::async, [&A,&b] { return solve(A,b); });
      
      while (1) {
	if(kbhit() && getchar()==27){ printf("\b\bcontsh");mkcont(argc,argv);system("bash");contop(argc,argv);system("rm .cont");}
	auto result = f.wait_for(chrono::seconds(1));
	if ( result != future_status::timeout) break;
      }
      x = f.get();
    }
    for(i=1;i<=m;i++){ Ux[i] = x[i]; Uy[i] = x[i+m];}

    if (defop("-o")) setanimefilename(getop("-o").c_str());
    plotncpolynomial1(Mid, x);
  }
  return 0;
}
