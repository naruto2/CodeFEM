#ifndef _EST_NAVIERSTOKES_HPP_
#define _EST_NAVIERSTOKES_HPP_
int navierstokes_init(const char *filename, double Re, double dt,
		      double (*u)(double x, double y),
		      double (*v)(double x, double y));
void navierstokes(sparse::matrix<double>&A, vector<double>&U,
		  vector<double>&b);
void plotuv(vector<double>&U);
void fprintuv(vector<double>&U);

int islabel(const char *label);
#endif
