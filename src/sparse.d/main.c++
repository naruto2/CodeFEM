#include "sparse.hpp"

using namespace sparse;

int main(){
  matrix<double> A(3);

  A[1][1] = 1.1;
  A[2][2] = 2.2;
  printmatrix(A);
  return 0;
}
