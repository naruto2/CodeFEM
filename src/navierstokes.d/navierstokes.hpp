#ifndef _EST_NAVIERSTOKES_HPP_
#define _EST_NAVIERSTOKES_HPP_
void navierstokes_init(const char *filename, double Re, double dt);
void navierstokes(sparse::matrix<double>&A, vector<double>&U, vector<double>&b);
void plotuv(vector<double>&U);
void fprintuv(vector<double>&U);
#endif
