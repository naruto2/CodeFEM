#include <cstdlib>
#include "solver.hpp"
#include "cg.h"
#include "bicgstab.h"
#include "cgs.h"
#include "bicg.h"
#include "qmr.h"
#include "jacobi.h"
#include "incholesky.h"
#include "Preconditioner.hpp"
#include "psc98.hpp"
#include "est/matrix.hpp"

void blockmatrix(matrix &A, matrix &B);

/*!
 * LU分解(ピボット交換なし)
 *  - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分解する
 *  - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
 *  - LとUを一つの行列にまとめた形で結果を返す
 * @param[inout] A n×nの係数行列．LU分解した結果を格納する．
 * @param[in] n 行列の大きさ
 * @return 1:成功,0:失敗
 */
int LUDecomp(matrix &A, int n)
{
  if(n <= 0) return 0;

  for(int i = 0; i < n; ++i){
    // l_ijの計算(i >= j)
    for(int j = 0; j <= i; ++j){
      double lu = A[i][j];
      for(int k = 0; k < j; ++k){
	lu -= A[i][k]*A[k][j];    // l_ik * u_kj
      }
      A[i][j] = lu;
    }
    
    // u_ijの計算(i < j)
    for(int j = i+1; j < n; ++j){
      double lu = A[i][j];
 
      for(int k = 0; k < i; ++k){
	lu -= A[i][k]*A[k][j];    // l_ik * u_kj
      }
      A[i][j] = lu/A[i][i];
    }
  }

  return 1;
}

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
  matrix A;
  vector<double> b;
  getprob(A,b);
  int n = A.size();
  vector<double> x(n);
  Preconditioner M, M2;
  int max_iter = 100000;
  double tol = 0.000001;

  M.ic(A);
  //BiCGSTAB(A, x, b, M, max_iter, tol);
  A.T(); QMR(A, x, b, M, M2, max_iter, tol);
  check(x);
  return 0;
}
