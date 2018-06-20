#include <vector>
#include "est/sparse.hpp"
#include "est/psc98.hpp"

extern "C" {
  int main2(void);
};


int main(int argc, char **argv){

  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  /* x = main2() */
  main2();
  psc98_check(x);
  return 0;
}
