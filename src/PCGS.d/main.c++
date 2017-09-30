#include <cstdio>
#include <cstdlib>
#include "est/sparse.hpp"
#include <f2c.h>

extern "C" {
  int pcgs_(doublereal *d__, doublereal *a, integer *ia,
	    integer *n, integer *n1, integer *nl,
	    doublereal *b, doublereal *eps,
	    integer *itr, doublereal *s, doublereal *x, doublereal *dd,
	    doublereal *p, doublereal *q, doublereal *r__, doublereal *r0,
	    doublereal *e, doublereal *h__, doublereal *w, integer *m,
	    integer *ier);
};

extern "C" {
  void  estiva_ary2(void**, long, long, size_t);
}

int PCGS(sparse::matrix<double>&matrix_A, vector<double>&vector_b)
{
  doublereal *D; doublereal   *A;  integer   *IA;
  integer     N; integer      N1;  integer    NL;
  doublereal *B; doublereal  EPS;
  integer   ITR; doublereal    S;  doublereal *X; doublereal *DD;
  doublereal *P; doublereal   *Q;  doublereal *R; doublereal *R0;
  doublereal *E; doublereal   *H;  doublereal *W; integer     *M;
  integer   IER;

  N  = matrix_A.size()-1;
  printf("N=%d\n",N);

  EPS = 0.001;
  ITR = 1000000;
  S   = 0.0;  
  int i, j, k, n;

  /* calc NL */
  NL=0;
  for (i=1; i<matrix_A.size(); i++){
    n = 0;
    for ( auto it: matrix_A[i] ) {
      j = it.first;
      if ( i < j && matrix_A[i][j] != 0.0 ) n++;
    }
    NL = max(NL,n);
  }  
  for (i=1; i<matrix_A.size(); i++){
    n = 0;
    for ( auto it: matrix_A[i] ) {
      j = it.first;
      if ( i > j && matrix_A[i][j] != 0.0 ) n++;
    }
    NL = max(NL,n);
  }  
  printf("NL=%d\n",NL);
  N1 = N+2*NL;

  D  = (doublereal*)calloc(sizeof(doublereal),N); // 一次元配列 D(N)
  B  = (doublereal*)calloc(sizeof(doublereal),N); // 一次元配列 B(N)
  A  = (doublereal*)calloc(sizeof(doublereal),N1*2*NL);// 二次元配列 A(N1,2*NL)
  IA = (integer*)calloc(sizeof(integer),N1*2*NL);      // 二次元配列IA(N1,2*NL)
  R  = (doublereal*)calloc(sizeof(doublereal),N);      // 一次元配列 R(N)

  X  = (doublereal*)calloc(sizeof(doublereal),N+1); //一次元配列で要素は0-N
  DD = (doublereal*)calloc(sizeof(doublereal),N+1); //一次元配列で要素は0-N
  P  = (doublereal*)calloc(sizeof(doublereal),N+1); //一次元配列で要素は0-N
  Q  = (doublereal*)calloc(sizeof(doublereal),N+1); //一次元配列で要素は0-N

  M  = (integer*)calloc(sizeof(integer),2*N); // 作業用
  
  R0 = (doublereal*)calloc(sizeof(doublereal),N); // 要素数は N
  E  = (doublereal*)calloc(sizeof(doublereal),N); // 要素数は N
  H  = (doublereal*)calloc(sizeof(doublereal),N); // 要素数は N

  W  = (doublereal*)calloc(sizeof(doublereal),N+1); //　要素は 0-N


  /* set D */
  for(i=1; i<matrix_A.size(); i++) D[i-1] = matrix_A[i][i];

  static doublereal **estivaA;
  static integer    **estivaIA;
  estiva_ary2((void**)&estivaA,2*NL, N+2*NL,sizeof(doublereal));
  estiva_ary2((void**)&estivaIA,2*NL, N+2*NL,sizeof(integer));
  
  
  for(i=1; i<=N; i++)for(k=0, j=1; j<i; j++)if(matrix_A[i][j] != 0.0){
	estivaA[k][i-1]  = matrix_A[i][j]  ;
	estivaIA[k][i-1] = j       ;
	k++                  ;
      }

  for(i=1; i<=N; i++)for(k=NL, j=i+1; j<=N; j++)if(matrix_A[i][j] != 0.0){
	estivaA[k][i-1]  = matrix_A[i][j] ;
	estivaIA[k][i-1] = j       ;
	k++                  ;
      }
  
  /* set B */
  for ( i=1; i<vector_b.size(); i++) B[i-1] = vector_b[i];
  
  printf("B ");for (i=0;i<4;i++) printf("%f ",B[i]);
  printf("\n");
  pcgs_(D,estivaA[0],estivaIA[0],&N,&N1,&NL,B,&EPS,&ITR,&S,X,DD,P,Q,R,R0,E,H,W,M,&IER);
  printf("X ");for (i=0;i<4;i++) printf("%f ",X[i+1]);
  printf("\n");



  free(R0);
  free(E);
  free(H);
  free(W);
  free(D);
  free(DD);
  free(P);
  free(Q);
  free(R);
  free(M);
  free(A);
  free(IA);
  free(B);
  free(X);
}

int PCGS(sparse::matrix<double>&A, vector<double>&b);
int GSLV1(sparse::matrix<double>&A, vector<double>&b);

int main()
{
  sparse::matrix<double> A(5);
  vector<double> b(5);
  A[1][1] = 1.0;  A[1][2] = 1.0;  A[1][3] = 1.0;  A[1][4] = 1.0;
  A[2][1] = 1.0;  A[2][2] = 1.0;  A[2][3] = 1.0;  A[2][4] =-1.0;
  A[3][1] = 1.0;  A[3][2] = 1.0;  A[3][3] =-1.0;  A[3][4] = 1.0;
  A[4][1] = 1.0;  A[4][2] =-1.0;  A[4][3] = 1.0;  A[4][4] = 1.0;

  b[1] = 0.0;
  b[2] = 4.0;
  b[3] =-4.0;
  b[4] = 2.0;

  
  PCGS(A,b);
  return 0;
}
