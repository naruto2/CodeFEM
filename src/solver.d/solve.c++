#include <vector>
#include "est/solver.hpp"
#include "est/ViennaCL.hpp"
#include "est/op.hpp"

using namespace std;

vector<double> solve(matrix &A, vector<double> &b) {
  Preconditioner M;
  //M.ic(A);


  if ( getop("-solver") == "bicgstab")   { printf("bicgstab\n");    return bicgstab(M,A,b); }
  if ( getop("-solver") == "cgs")        { printf("cgs\n");         return cgs(M,A,b); }
  if ( getop("-solver") == "gpubicgstab"){ printf("gpubicgstab\n"); return gpubicgstab(A,b); }
  return gpubicgstab(A,b);
}
