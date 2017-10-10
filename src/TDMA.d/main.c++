#include <stdio.h>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>


vector<double> TDMA(sparse::matrix<double>&A,vector<double>&d);

  
int main(int argc, char **argv){
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
  return 0;
}
