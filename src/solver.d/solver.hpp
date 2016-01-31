#ifndef _EST_SOLVER_HPP_
#define _EST_SOLVER_HPP_

#include <cstdio>
#include <cmath>
#include "matrix.hpp"


double dot(const vector<double>&x, const vector<double>&y);
double norm(const vector<double>&x);
vector<double>& operator-(const vector<double>&x, vector<double>&y);
vector<double>& operator+(vector<double>&x, vector<double>&y);
vector<double>& operator*(double a, vector<double>&x);
vector<double>& operator*(const matrix&A, vector<double>&b);

typedef double Real;
typedef matrix Matrix;
typedef vector<double> Vector;

int IncompleteCholeskyDecomp2(matrix &A, matrix &L, vector<double> &d, int n);

class Preconditioner{
public:

  mutable matrix A, L;
  mutable vector<double> d;
  mutable int solver = 0;
  vector<double>& nsolve(vector<double>& p)const;
  void jacobi(matrix& pA);
  vector<double>& scsolve(vector<double>& p) const;
  vector<double>& icsolve(vector<double>& b) const;
  vector<double> & solve(vector<double>& b) const;
  vector<double>& trans_solve(vector<double>& p) const;
  void ic(matrix& Ap){
    solver = 1;
    A = Ap;
    int n = A.size();
    L.resize(n);
    d.resize(n);
    A.sync();
    IncompleteCholeskyDecomp2(A, L, d, n);
  }
};

vector<double> cg(Preconditioner& M, matrix& A, vector<double>& b);
vector<double> ir(Preconditioner& M, matrix& A, vector<double>& b);
vector<double> cgs(Preconditioner& M, matrix& A, vector<double>& b);
vector<double> bicgstab(Preconditioner& M, matrix& A, vector<double>& b);
vector<double> bicg(Preconditioner& M, matrix& A, vector<double>& b);
vector<double> qmr(Preconditioner& M, Preconditioner& M2, matrix& A, vector<double>& b);
vector<double> gmres(Preconditioner& M, matrix& A, vector<double>& b);
vector<double> cheby(Preconditioner& M, matrix& A, vector<double>& b, double mine, double maxe);
vector<double> jacobi(Preconditioner& M, matrix& A, vector<double>& b);

void blockmatrix(matrix &A, matrix &B);

#endif
