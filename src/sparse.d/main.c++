#include "sparse.hpp"

using namespace sparse;

int main(){
  matrix<double> A(11);

  A[1][1] = 1.1;
  A[2][2] = 2.2;
  A[9][9] = 9.9;
  plotmatrix(A);
  return 0;
}
