#include <cstdlib>
#include "est/sparse.hpp"
#include "est/solver.hpp"
#include "est/psc98.hpp"

int main(int argc, char **argv){
  
#if 0
  matrix B(3);
  B[0][0] = 1.00; B[0][1] = 0.50; B[0][2] = 0.00; 
  B[1][0] = 0.50; B[1][1] = 1.00; B[1][2] = 0.50; 
  B[2][0] = 0.00; B[2][1] = 0.50; B[2][2] = 1.00;
  LUDecomp(B,B.size());
  cout << B << endl;
  return 0;
#endif
  sparse::matrix<double> A;
  vector<double> b;
  psc98_init(A,b);
  Preconditioner M, M2;

  M.ic(A); 
  //  vector<double> x = cheby(M,A,b,19.738218,129105.765540);
  vector<double> x = cg(M,A,b);
  psc98_check(x);
  return 0;
}

