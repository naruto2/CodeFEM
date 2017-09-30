#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"
#include "estiva/eblas.h"


typedef struct {
  long   m, n, *col_ind, *row_ptr, *diag_ptr;
  double *val, *pivots;
} CRS;

typedef struct{
  double **A;
  long   **IA;
  long   n, w, I, J;
  double a;
} MX;



extern "C" {
  int   estiva_pcgssolver(void *Apointer, double *xk, double *b);
  void  estiva_ary1(void**, long, size_t);
  void  estiva_initmx(MX **Ap, long i, long j);
  double *estiva_mx(MX *T, long i, long j);
};

int main()
{
  static MX *A;
  static double *x, *b;


  estiva_ary1((void**)&x,5,sizeof(x));
  estiva_ary1((void**)&b,5,sizeof(b));
  estiva_initmx((MX**)&A,5,5);
  
  (*estiva_mx(A,1,1)) = 1.0;  (*estiva_mx(A,1,2)) = 1.0;  (*estiva_mx(A,1,3)) = 1.0;  (*estiva_mx(A,1,4)) = 1.0;
  (*estiva_mx(A,2,1)) = 1.0;  (*estiva_mx(A,2,2)) = 1.0;  (*estiva_mx(A,2,3)) = 1.0;  (*estiva_mx(A,2,4)) =-1.0;
  (*estiva_mx(A,3,1)) = 1.0;  (*estiva_mx(A,3,2)) = 1.0;  (*estiva_mx(A,3,3)) =-1.0;  (*estiva_mx(A,3,4)) = 1.0;
  (*estiva_mx(A,4,1)) = 1.0;  (*estiva_mx(A,4,2)) =-1.0;  (*estiva_mx(A,4,3)) = 1.0;  (*estiva_mx(A,4,4)) = 1.0;

  b[1] = 0.0;
  b[2] = 4.0;
  b[3] =-4.0;
  b[4] = 2.0;
  
  estiva_pcgssolver(A,x,b);
  for(int i=1;i<5; i++) printf("%f ",x[i]);
  return 0;
}
