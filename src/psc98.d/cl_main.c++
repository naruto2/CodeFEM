#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/psc98.hpp>

vector<double> TDMA(sparse::matrix<double>&A,vector<double>&d);

int isTridiagonal(sparse::matrix<double>&A)
{
  int n = A.size();
  for (int i=1; i<n; i++)
    for ( auto it : A[i] ){ int j = it.first;
      if ( abs(i-j)>1 ) return 0;
    }
  return 1;
}

int main(int argc, char **argv){
  cl_bicgstab_init(argc,argv);

  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  if (isTridiagonal(A))
    x = TDMA(A,b);
  else
    x = cl_bicgstab(A,b);
  psc98_check(x);
  return 0;
}
