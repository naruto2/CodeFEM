#include <vector>
#include "est/sparse.hpp"
#include "est/psc98.hpp"

extern int PCGS(sparse::matrix<double>&AA, vector<double>&xx, vector<double>&bb);
  
int main(int argc, char **argv){

  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  x.resize(A.size());
  PCGS(A,x,b);
  psc98_check(x);
  return 0;
}
