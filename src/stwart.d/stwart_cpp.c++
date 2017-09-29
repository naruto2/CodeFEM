#include <cstdlib>
#include "est/sparse.hpp"
#include <f2c.h>

extern "C" {
  int acutualmain(int argc, char **argv, integer *ib, integer l, integer *mb, integer m);

  int stwart_(integer *, integer *, integer *, integer *, integer *,
	      integer *, integer *, integer *, integer *, integer *,
	      integer *, integer *, integer *, integer *, integer *,
	      integer *, integer *);
};


integer *generate_ma(sparse::matrix<double>&A)
{
  int i, j, n;
  integer *ma;
  ma = (integer*)calloc(sizeof(integer),A.size());
  
  n = 0;
  for ( i = 1; i < A.size(); i++ ) {
    for ( auto it:A[i]) {
      j = it.first;
      if ( A[i][j] != 0.0 ) n++;
    }
    ma[i] = n;
  }
  return ma;
}

integer generate_ia(sparse::matrix<double>&A, integer **ia)
{
  int i, j, n;

  n = 0;
  for ( i = 1; i < A.size(); i++ ) {
    for ( auto it:A[i] ) {
      j = it.first;
      if ( A[i][j] != 0.0 ) n++;
    }
  }
  *ia = (integer*)calloc(sizeof(integer),n);

  n = 0;
  for ( i = 1; i < A.size(); i++ ) {
    for ( auto it:A[i] ) {
      j = it.first;
      if ( A[i][j] != 0.0 ) (*ia)[n++] = j;
    }
  }
  return n;
}


int stwart(sparse::matrix<double>&A)
{
    integer *ma, m;
    ma = generate_ma(A);
    m  = A.size()-1; 
    
    integer *ia, l;
    l = generate_ia(A,&ia);
    
      integer *R    = (integer*)calloc(sizeof(integer),m);
      integer *C    = (integer*)calloc(sizeof(integer),m);
      integer *IR   = (integer*)calloc(sizeof(integer),m);
      integer *IC   = (integer*)calloc(sizeof(integer),m);
      integer *JROW = (integer*)calloc(sizeof(integer),m);
      integer *JCOL = (integer*)calloc(sizeof(integer),m);
      integer *IP   = (integer*)calloc(sizeof(integer),m);
      integer *JP   = (integer*)calloc(sizeof(integer),m);
      integer *IW   = (integer*)calloc(sizeof(integer),m);

      integer kerns=0, mend=0, lg=0, ier=0;
      
      printf("stwart\n");
      stwart_(ia, &l, ma, &m, R, C, IR, IC, JROW, JCOL, IP, JP, &kerns,
                &mend, IW, &lg, &ier);

      printf("m=%d, l=%d, kerns=%d, mend=%d, lg=%d, ier=%d\n",
              m,    l,    kerns,    mend,    lg,    ier);

      printf("R: ");
      for(int i=0; i<m; i++) printf("%d ",R[i]);
      printf("\n");

      printf("C: ");
      for(int i=0; i<m; i++) printf("%d ",C[i]);
      printf("\n");

      printf("IW: ");
      for(int i=0; i<m; i++) printf("%d ",IW[i]);
      printf("\n");

      printf("JCOL: ");
      for(int i=0; i<m; i++) printf("%d ",JCOL[i]);
      printf("\n");

      printf("JROW: ");
      for(int i=0; i<m; i++) printf("%d ",JROW[i]);
      printf("\n");

      free(R);
      free(C);
      free(IR);
      free(IC);
      free(JROW);
      free(JCOL);
      free(IP);
      free(JP);
}
