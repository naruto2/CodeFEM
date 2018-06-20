#include <stdio.h>
#include "est/sparse.hpp"
#include <f2c.h>


extern "C" {
  int pcgs_(doublereal *, doublereal *, integer *,
	    integer *, integer *, integer *, doublereal *, doublereal *,
	    integer *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
};


static int memory_error(void)
{
  fprintf(stderr,"PCGS() can't malloc()\n");
  return 1;
}


int PCGS(sparse::matrix<double>&AA, vector<double>&xx, vector<double>&bb)
{
  /* (i)  引数の型と種類                                                */
  doublereal*  D, *B      ;/*    D  一次元配列 D(N), B(N)                */
  doublereal*  A         ;/*    D  二次元配列 A(N1, 2 * NL)             */
  integer*    IA         ;/*    I  二次元配列 IA(N1, 2 * NL)            */
  doublereal*  R         ;/*    D  一次元配列 R(N)                      */
  integer     NL, N1, N, ITR, IER                                       ;
  doublereal EPS, S                                                     ;
  doublereal *X,*DD,*P,*Q;/*    D  一次元配列で, 要素は 0〜N            */
  integer*    M          ;/*    I  一次元配列  M(2 * N)  作業用         */

  /*--    DからRまでと, EPSからSまでの引数は, サブルーチンPCGと同じ   --*/
  doublereal *R0, *E,* H ;/*    D  一次元配列  要素数はN                */
  doublereal *W          ;/*    D  一次元配列  W(0:N)                   */

  integer   i, j, k      ;

  N  = (integer)AA.size()-1;

  for(NL=0, i=1; i<=N; i++){
    k = 0;
    for ( auto AAi: AA[i] ) k++;
    if(NL<k) NL = k;
  }

  B = (doublereal*)calloc(sizeof(doublereal),N); if (!B) return memory_error();
  D = (doublereal*)calloc(sizeof(doublereal),N); if (!D) return memory_error();
  A = (doublereal*)calloc(sizeof(doublereal),2*NL*(N+2*NL));
  if (!A) return memory_error();
  IA= (integer*   )calloc(sizeof(integer   ),2*NL*(N+2*NL));
  if (!IA) return memory_error();
  R = (doublereal*)calloc(sizeof(doublereal),N); if (!R) return memory_error();

  X = (doublereal*)calloc(sizeof(doublereal),N+1);if(!X) return memory_error();
  DD= (doublereal*)calloc(sizeof(doublereal),N+1);if(!DD)return memory_error();
  P = (doublereal*)calloc(sizeof(doublereal),N+1);if(!P) return memory_error();
  Q = (doublereal*)calloc(sizeof(doublereal),N+1);if(!Q) return memory_error();

  M = (integer*   )calloc(sizeof(integer   ),2*N);if(!M) return memory_error();

  R0= (doublereal*)calloc(sizeof(doublereal),N);if(!R0)return memory_error();
  E = (doublereal*)calloc(sizeof(doublereal),N);if(!E) return memory_error();
  H = (doublereal*)calloc(sizeof(doublereal),N);if(!H) return memory_error();

  W = (doublereal*)calloc(sizeof(doublereal),N+1);if(!W) return memory_error();
  

  /* (ii) 主プログラム → サブルーチン                                  */
  /*      サブルーチンをCALLするときには, つぎの値を与える.             */
  /*  D   : 配列Dの第1〜第n位置に行列Aの対角要素を入れておく            */
  /*  N   : 行列Aの行数を入れておく.                                    */
  /*  N1  : 配列Aの行数を入れておく. N1≧N+2*NL でないといけない.       */
  /*  NL  : 行列Aの各行における非ゼロ要素数の最大値を入れておく.        */
  /*  B   : 連立一次方程式の右辺を入れておく.                           */
  /* EPS  : 収束判定置を入れておく. ふつうは 1.×10^(-7)                */
  /* ITR  : 打切りまでの最大繰返し回数を入れておく.                     */

  /*--    D, N, N1, NL, B, EPS, ITR については, サブルーチンPCGと同じ --*/

  /*  A   : 配列Aの各行1〜NL要素は, 行列Aの下三角部分の各行の非ゼロ     */
  /*        要素を入れる. また, 各行のNL+1〜2*NL要素は, 上三角部分      */
  /*        の各行の非ゼロ要素を入れる. ただし, 各行の対角要素は配列Dに */
  /*        入れる.                                                     */
  /* IA   : 配列Aに入れた要素の列番号を, 対応する位置に入れておく.      */
  /*  S   : LUCGS法のとき 0., MLUCG法のときσ(>0)を入れておく.          */

  for(i=1; i<=N; i++) D[i-1] = AA[i][i];

  for(i=1; i<=N; i++) {
    k = 0;
    for ( auto AAi : AA[i] ) if ( AAi.first<i ) {
	j = AAi.first;
	A[k*2*NL+(i-1)]  = AA[i][j];
	IA[k*2*NL+(i-1)] = j;
	k++;
      }
  }

  
  for(i=1; i<=N; i++) {
    k = NL;
    for ( auto AAi : AA[i] ) if ( i<AAi.first ) {
	j = AAi.first;
	A[k*2*NL+(i-1)]  = AA[i][j];
	IA[k*2*NL+(i-1)] = j;
	k++;
      }
  }

  N1  =  N+2*NL;

  for(i=1; i<=N; i++) B[i-1] = bb[i];
  EPS = 1.0e-7 ;
  ITR = N      ;
  S   = 0.     ;

  /* (a) リンクの方法                                                   */
  /* CALL PCGS(D,A,IA,N,N1,NL,B,EPS,ITR,S,X,DD,P,Q,R,R0,E,H,W,M,IER)    */

  pcgs_(D,A,IA,&N,&N1,&NL,B,&EPS,&ITR,
       &S,X,DD,P,Q,R,R0,E,H,W,M,&IER);
  
  for(i=1; i<=N; i++) xx[i] = X[i-1];
  printf("ITR = %ld\n",ITR);

  /*  free(B);
  free(D);
  free(A);
  free(IA);
  free(R);
  free(X);
  free(DD);
  free(P);
  free(Q);
  free(M);
  free(R0);
  free(E);
  free(H);
  free(W);

  */  return IER;
}
