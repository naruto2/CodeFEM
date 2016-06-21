#include <vector>
#include "est/solver.hpp"
#include "est/ViennaCL.hpp"

using namespace std;

vector<double> solve(matrix &A, vector<double> &b) {
  Preconditioner M;
  //M.ic(A);

  return bicgstab(M,A,b);
  return gpusolver(A,b);
  return cgs(M,A,b);
}
