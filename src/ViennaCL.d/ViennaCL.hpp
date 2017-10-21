#ifndef _EST_VIENNACL_HPP_
#define _EST_VIENNACL_HPP_
#define VIENNACL_WITH_OPENCL

vector<double> vcl_cg(sparse::matrix<double>& A, vector<double>& b);
vector<double> vcl_cgilut(sparse::matrix<double>& A, vector<double>& b);
vector<double> vcl_bicgstab(sparse::matrix<double>& A, vector<double>& b);
vector<double> vcl_gmres(sparse::matrix<double>& A, vector<double>& b);

#endif
