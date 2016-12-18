// Ax=b を解く

#include <iostream>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;

#include "EigenTools.h"
#include "bicgstab.h"


// Smatrix x = Solve(A,b); LU分解による連立一次方程式の求解
Vector Solve(Smatrix &A, Vector &b) {

  SparseLU<Smatrix> Solver;
  Solver.compute(A);
  if(Solver.info() != Success) {
    Error("decomposition failed");
    Vector z;
    return z;
  }

  Vector x = Solver.solve(b);
  if(Solver.info() != Success) Error("solving failed");

  return x;
}


// Smatrix x = Isolve(A,b); BiCGSTAB法による連立一次方程式の求解
Vector Isolve(Smatrix &A, Vector &b) {

  BiCGSTAB<Smatrix> Solver;
  Solver.compute(A);
  if(Solver.info() != Success) {
    Error("decomposition failed");
    Vector z;
    return z;
  }

  Vector x = Solver.solve(b);
  if(Solver.info() != Success) Error("solving failed");

  cout << "#iterations:     " << Solver.iterations() << endl;
  cout << "estimated error: " << Solver.error()      << endl;
  
  return x;
}

Vector Ssolve(Smatrix &A, Vector &b) {

  printf("size = %d\n",b.size());
  printf("row = %d\n",A.rows());
  printf("col = %d\n",A.cols());

  Smatrix AT = A.transpose();
  Smatrix A11 = A + AT; Smatrix A12 = A - AT;
  Smatrix A21 = AT -A ; Smatrix A22 = -A - AT;

  Smatrix C, D, E;

  int n = A.cols();

  C.resize(n, 2*n);
  C.middleCols(0,A11.cols()) = A11;
  C.middleCols(A11.cols(),A12.cols()) = A12;

  D.resize(n, 2*n);
  D.middleCols(0,A21.cols()) = A21;
  D.middleCols(A21.cols(),A22.cols()) = A22;

  Smatrix CT = C.transpose();
  Smatrix DT = D.transpose();

  E.resize(2*n, 2*n);
  E.middleCols(0,CT.cols()) = CT;
  E.middleCols(CT.cols(),DT.cols()) = DT;
  
  cout << E << endl;

  Vector c = -b;
  Vector bb(b.rows()+c.rows());
  
  bb << 2*b, 2*c;

  cout << bb << endl;

  return Solve(E,bb);
}

// メイン関数
int main(void){

  const double s=19, u=21, p=16, e=5, r=18, l=12;

  Tri(0,0, s);
  Tri(1,0, l);
  Tri(4,0, l);
  Tri(1,1, u);
  Tri(2,1, l);
  Tri(4,1, l);
  Tri(0,2, u);
  Tri(2,2, p);
  Tri(0,3, u);
  Tri(3,3, e);
  Tri(3,4, u);
  Tri(4,4, r);

  Smatrix A = MapSmatrix(5,5);

  double  B[] = {1,1,1,1,1}; // 右辺値（普通のC++配列で与える例）
  Vector   b  = MapVector(B, 5);

  Vector x = Isolve(A,b);

  cout << "solution Vector" << endl << x << endl;

  return 0;
}
