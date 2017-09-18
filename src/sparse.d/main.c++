#include "sparse.hpp"

using namespace sparse;

int main(){
  matrix<double> A;

  A.resize(3);
  A[1][1] = 10.1;
  A[2][2] = 20.1;
  printmatrix(A);
  return 0;
}
