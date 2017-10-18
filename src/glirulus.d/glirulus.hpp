#ifndef _EST_GLIRULUS_HPP_
#define _EST_GLIRULUS_HPP_

#include <cstdio>
#include <cmath>
#include "mmio.h"
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/op.hpp>
#include <vector>
#include "est/GLU1.hpp"

double L_inf(vector<double>&x);
int glirulus_init(int argc, char **argv);
int glirulus_mm(sparse::matrix<double>&A,vector<double>&b);
vector<double> glirulus(sparse::matrix<double>&A,vector<double>&b);
int glirulus_check(sparse::matrix<double>&A,vector<double>&x,vector<double>&b);

#endif
