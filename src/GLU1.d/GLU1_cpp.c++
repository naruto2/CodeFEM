#include <stdio.h>
#include "est/sparse.hpp"
#include <f2c.h>

extern "C" {
  int glu1_(doublereal *, integer *, integer *,
	    doublereal *, doublereal *, integer *, integer *);
};
  

int GLU1(sparse::matrix<double>&A)
{
  integer N=A.size()-1, IER=0;
  doublereal EPS = -1.0;
  
  integer    *IP = (integer*)calloc(sizeof(integer),N);
  doublereal *WK = (doublereal*)calloc(sizeof(doublereal),N);
  doublereal *a  = (doublereal*)calloc(sizeof(doublereal),N*N);

  for (int i = 1; i < A.size(); i++ ) for ( auto it:A[i] ) {
      int j = it.first;
      a[ (i-1) + N*(j-1) ] = A[i][j];
    }

  glu1_(a, &N, &N, &EPS, WK, IP, &IER);

  free(IP);
  free(WK);
  free(a);
  return (int)IER;
}
