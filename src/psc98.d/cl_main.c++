#include <cstdlib>
#include "est/sparse.hpp"

void getprob(sparse::matrix<double>&A, vector<double>& b);
void check(vector<double>& x);
vector<double> sparse__bicgstab(sparse::matrix<double>& A, vector<double>& b);
void cl_init(int argc, char **argv);

int main(int argc, char **argv){
  cl_init(argc,argv);
  sparse::matrix<double> A;
  vector<double> b;
  getprob(A,b);
  vector<double> x;

  int max_iter = 100000;
  double tol = 0.000001;

  x = sparse__bicgstab(A,b);

  check(x);
  return 0;
}
