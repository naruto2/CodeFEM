#ifndef _EST_BICGSTAB_HPP_
#define _EST_BICGSTAB_HPP_

int cl_bicgstab_init(int argc,char **argv);
vector<double> cl_bicgstab(sparse::matrix<double>& A, vector<double>& b);
  
#endif
