#include <stdio.h>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>

vector<double> TDMA(sparse::matrix<double>&A,vector<double>&d)
{
  int i, n = A.size();

  vector<double> a(n), b(n), c(n), P(n), Q(n), x(n);
  
  for ( i=1; i<n; i++)   a[i] =  A[i][i];
  for ( i=1; i<n-1; i++) b[i] = -A[i][i+1];
  for ( i=2; i<n;   i++) c[i] = -A[i][i-1];

  P[1] = b[1]/a[1];
  Q[1] = d[1]/a[1];

  for( i=2;i<n; i++) P[i] = b[i]/(a[i]-c[i]*P[i-1]);
  for( i=2;i<n; i++) Q[i] = (d[i]+c[i]*Q[i-1])/(a[i]-c[i]*P[i-1]);

  x[n-1] = Q[n-1];
  for ( i=n-2; i>0; i--) x[i] = P[i]*x[i+1]+Q[i];

  return x;
}


main(int argc, char **argv){
  cl_bicgstab_init(argc,argv);
  int i, n = 512;

  sparse::matrix<double> A(n+1);
  vector<double> x(n+1), b(n+1);

  for ( i =1; i<=n; i++){
    A[i][i] = 2.0; b[i] = 1.0;
  }
  for (i=1;i<n;i++) A[i][i+1] = -1.0;
  for (i=1;i<n;i++) A[i+1][i] = -1.0;

  x = TDMA(A,b);

  for (i=1;i<=n;i++) printf("%f\n", x[i] );
}
