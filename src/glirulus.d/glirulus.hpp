#ifndef _EST_GLIRULUS_HPP_
#define _EST_GLIRULUS_HPP_

#include <cstdio>
#include <cmath>
#include <vector>

double L_inf(vector<double>&x);
int glirulus_init(int argc, char **argv);
int glirulus_mm(sparse::matrix<double>&A,vector<double>&b);
vector<double> glirulus(sparse::matrix<double>&A,vector<double>&b);
double glirulus_check(sparse::matrix<double>&A,vector<double>&x,const vector<double>&b);

#endif
