#include <vector>
#include "est/sparse.hpp"
#include "est/solver.hpp"
#include "est/ViennaCL.hpp"
#include "est/op.hpp"

using namespace std;

vector<double> solve(sparse::matrix<double> &A, vector<double> &b) {
  Preconditioner M;
  //M.ic(A);

#if 0
  if ( getop("-solver") == "bicgstab")   { printf("bicgstab\n");    return bicgstab(M,A,b); }
  if ( getop("-solver") == "cgs")        { printf("cgs\n");         return cgs(M,A,b); }
  if ( getop("-solver") == "vcl_bicgstab"){ printf("vcl_bicgstab\n"); return vcl_bicgstab(A,b); }
#endif
  return vcl_bicgstab(A,b);
}
