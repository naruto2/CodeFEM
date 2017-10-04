#include <cstdlib>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>

void getprob(sparse::matrix<double>&A, vector<double>& b);
void check(vector<double>& x);

int main(int argc, char **argv){
  cl_bicgstab_init(argc,argv);
  sparse::matrix<double> A;
  vector<double> b;
  getprob(A,b);
  vector<double> x;

  int max_iter = 100000;
  double tol = 0.000001;

  x = cl_bicgstab(A,b);

  check(x);
  return 0;
}
