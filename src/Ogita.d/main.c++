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

vector<xyc> ncpolynomial1(vector<xyc> Z, vector<nde> &N );
void estiva_forgammap1(long *x);
int estiva_forgammap1_loop(long *x, const char *name, vector<xyc> &Z);
void printmatrix(matrix &A, const char *name);  
void printvector(vector<double> &b, const char *name);
void plotncpolynomial1(vector<xyc> Mid, vector<double> x);
void setanimefilename(const char *fname);
void squaremesh(int n, vector<xyc> &Z);

void start_baton(void);
int end_baton(void);
#define batonth() for(start_baton();end_baton();)

void mkcont(int argc, char**argv);
void contop(int &argc0, char** &argv0);
int kbhit(void);

vector<double> S_(vector<xyc> &Z, vector<nde> &N);
matrix M__(vector<xyc> &Mid, vector<nde> &N, vector<double> &S);
matrix K__(vector<xyc> &Mid, vector<xyc> &Z, vector<nde> &N, vector<double> &S);
matrix Hx__(vector<xyc> &Mid, vector<xyc> &Z, vector<nde> &N);
matrix Hy__(vector<xyc> &Mid, vector<xyc> &Z, vector<nde> &N);
void A__(matrix &A, vector<xyc> &Mid, vector<nde> &N, matrix &M, double tau, matrix &K, matrix &Hx, matrix &Hy);
void Rhs(vector<double> &b, vector<xyc> &Mid, vector<nde> &N,matrix &M,double t,vector<double> &Fx,vector<double> &Fy,
	 vector<double> &Ux, vector<double> &Uy, vector<double> &x);
void boundary_condition(vector<nde> &N, vector<xyc> &Mid, matrix &A, vector<double> &b);

  
	 

int main(int argc, char ** argv)
{
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  {
    unsigned long n = 16;
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

  for ( k = 1; k<=60; k++ ) {
    A__(A, Mid, N, M,t,K,Hx,Hy);
    Rhs(b, Mid, N, M, t, Fx, Fy, Ux, Uy, x);
    boundary_condition(N,Mid,A,b);

    printf("\nk = %ld  \n",k);

    //x = solve(A,b);

#define nonzeroloop(A) for(i=0; i<A.size(); i++) for ( auto it : A[i] ) if( A[i][(j=it.first)] != 0.0 )
    
    { int i, j; double t;
      
      /*
      NUM=3;
      A.clear();
      A.resize(NUM+1);
      A[1][1] = 1.0; A[1][2] = 1.0; A[1][3] = 0;
      A[2][1] = 0.0; A[2][2] = 1.0; A[2][3] = 1;
      A[3][1] = 0.0; A[3][2] = 0.0; A[3][3] = 1;
      b[1] = 2;      b[2] = 1;
      */


      matrix B;
      B.clear();
      B.resize(2*NUM+1);

      for(i=1;i<=NUM;i++)for(j=1;j<=NUM;j++)if((t=A[i][j]+A[j][i])!=0.0) {
	    B[i][j] = t;
	    B[i+NUM][j+NUM] = -t;
	  }
      for(i=1;i<=NUM;i++)for(j=1;j<=NUM;j++)if((t=A[i][j]-A[j][i])!=0.0) {
	    B[i][j+NUM] = t;
	    B[i+NUM][j] = -t;
	  }
      B[0][0] = 1;

      vector<double> xx(2*NUM+1), bb(2*NUM+1); 
      for(i=1;i<=NUM;i++){ bb[i] = 2.0*b[i]; bb[i+NUM] = -2.0*b[i];}

      //printmatrix(B,"B");
      xx = solve(B,bb);
      exit(0);
    }
    
    
    for(i=1;i<m;i++){ Ux[i] = x[i]; Uy[i] = x[i+m];}

    
    if (defop("-o")) setanimefilename(getop("-o").c_str());
    plotncpolynomial1(Mid, x);
  }
  return 0;
}
