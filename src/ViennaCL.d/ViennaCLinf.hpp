#ifndef _EST_VIENNACL_HPP_
#define _EST_VIENNACL_HPP_
#define VIENNACL_WITH_OPENCL

#include <viennacl/linalg/jacobi_precond.hpp>
#include "cginf.hpp"
#include "bicgstabinf.hpp"
#include "gmresinf.hpp"
#include "togpu.hpp"

vector<double> gpucg(matrix& A, vector<double>& b)
{
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n,0);
  viennacl::vector<double>     bgpu(n), xgpu(n);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::jacobi_precond< gpumatrix >  vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
  viennacl::linalg::cg_tag  custom_cg(1e-5,1000000);

  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_cg, vcl_jacobi);

  gpuvector2vector(xgpu,x);

  return x;
}


vector<double> gpubicgstab(matrix& A, vector<double>& b)
{
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n,0);
  viennacl::vector<double>     bgpu(n), xgpu(n);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::jacobi_precond< gpumatrix >  vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
  viennacl::linalg::bicgstab_tag  custom_bicgstab(1e-5,1000000);

  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab, vcl_jacobi);

  gpuvector2vector(xgpu,x);

  return x;
}


vector<double> gpugmres(matrix& A, vector<double>& b)
{
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n,0);
  viennacl::vector<double>     bgpu(n), xgpu(n);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::jacobi_precond< gpumatrix >  vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
  viennacl::linalg::gmres_tag  custom_gmres(1e-5,10000000,50);

  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_gmres, vcl_jacobi);

  gpuvector2vector(xgpu,x);

  return x;
}

#endif
