#include <cstdio>
#include <vector>
#include "est/sparse.hpp"

extern "C"{
  int estiva_bicgstab(long *Ai, long *Aj,double *Aa, double *X, double *B);
};


int estiva__bicgstab(sparse::matrix<double>&A,vector<double>&bp)
{
  int n=A.size(), i, j, k;
 
  for (k=0, i=1; i<A.size(); i++) for (auto it : A[i])
				    if (A[i][it.first] != 0.0 ) k++;
  vector<long>   Ai(k+1),Aj(k+1);
  vector<double> Aa(k+1),X(n), B(n);

  for (k=0,i=1; i<A.size(); i++) for (auto it : A[i]){
      j = it.first;
      if (A[i][j] != 0.0 ) {
	Ai[k] = i; Aj[k] = j; Aa[k] = A[i][j];
	k++;
      }
    }
  Ai[k] = 0;
  
  for(i=1; i<n; i++) B[i]=bp[i];

  estiva_bicgstab(&Ai[0],&Aj[0],&Aa[0],&X[0],&B[0]);

  for(i=1; i<n; i++) bp[i]=X[i];
  return 0;
}


#if 0
int main()
{
  sparse::matrix<double> A(3);
  vector<double>  b(3);
  A[1][1] = 1.0; A[1][2] = 1.0;
  A[2][1] = 0.0; A[2][2] = 1.0;

  b[1] = 1.0;
  b[2] = 1.0;
  estiva__bicgstab(A,b);
  printf("%f\n",b[1]);
  printf("%f\n",b[2]);
  return 0;
}
#endif
