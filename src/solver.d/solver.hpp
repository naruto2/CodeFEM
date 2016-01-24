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

#endif
