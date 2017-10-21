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

int IncompleteCholeskyDecomp2(sparse::matrix<double> &A, sparse::matrix<double> &L, vector<double> &d, int n);

class Preconditioner{
public:

  mutable sparse::matrix<double> A, L;
  mutable vector<double> d;
  mutable int solver = 0;
  vector<double>& nsolve(vector<double>& p)const;
  void jacobi(matrix& pA);
  vector<double>& scsolve(vector<double>& p) const;
  vector<double>& icsolve(vector<double>& b) const;
  vector<double> & solve(vector<double>& b) const;
  vector<double>& trans_solve(vector<double>& p) const;
  void ic(sparse::matrix<double>& Ap){
    solver = 1;
    A = Ap;
    int n = A.size();
    L.resize(n);
    d.resize(n);

    IncompleteCholeskyDecomp2(A, L, d, n);
  }
};

vector<double> cg(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);
vector<double> ir(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);
vector<double> cgs(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);
vector<double> bicgstab(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);
vector<double> bicg(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);
vector<double> qmr(Preconditioner& M, Preconditioner& M2, sparse::matrix<double>& A, vector<double>& b);
vector<double> gmres(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);
vector<double> cheby(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b, double mine, double maxe);
vector<double> jacobi(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b);

void blockmatrix(sparse::matrix<double> &A, sparse::matrix<double> &B);

vector<double> solve(sparse::matrix<double> &A, vector<double> &b);

#endif

