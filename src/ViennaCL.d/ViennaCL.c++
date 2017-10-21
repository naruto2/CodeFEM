#define VIENNACL_WITH_OPENCL

#include <vector>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <viennacl/linalg/ichol.hpp>
#include <viennacl/linalg/detail/ilu/chow_patel_ilu.hpp>
#include "est/sparse.hpp"
#include "cginf.hpp"
#include "bicgstabinf.hpp"
#include "gmresinf.hpp"
#include "togpu.hpp"
#include "est/op.hpp"

vector<double> vcl_cg(sparse::matrix<double>& A, vector<double>& b)
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


vector<double> vcl_cgilut(sparse::matrix<double>& A, vector<double>& b)
{
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n,0);
  viennacl::vector<double>     bgpu(n), xgpu(n);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);
  viennacl::linalg::cg_tag  custom_cg(1e-5,1000000);
  
  // viennacl::linalg::ichol0_tag ichol0_conf;
  // typedef viennacl::linalg::ichol0_precond<
  // viennacl::compressed_matrix<double> >  vcl_ilut_t;
  // vcl_ilut_t vcl_ilut(Agpu,ichol0_conf);
  // viennacl::linalg::jacobi_precond< gpumatrix >  vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());

  viennacl::linalg::chow_patel_tag icc_conf;
  typedef viennacl::linalg::chow_patel_icc_precond<
    viennacl::compressed_matrix<double> > vcl_icc_t;
  vcl_icc_t vcl_icc(Agpu,icc_conf);
  
  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_cg,vcl_icc);

  gpuvector2vector(xgpu,x);

  return x;
}


vector<double> vcl_bicgstab(sparse::matrix<double>& A, vector<double>& b)
{
  A[0][0] = 1.0; b[0] = 0.0;
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n,0);
  viennacl::vector<double>     bgpu(n), xgpu(n);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::bicgstab_tag  custom_bicgstab(1e-7,10000000);

  if ( getop("-pre") == "jacobi") {
    viennacl::linalg::jacobi_precond< gpumatrix >
      vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
  
    xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab, vcl_jacobi);
  }
  if ( getop("-pre") == "ilut" ) {
    viennacl::linalg::ilut_tag ilut_conf(10, 1e-7);
    //10 entries, rel. tol. 1e-7

    viennacl::linalg::ilut_precond< viennacl::compressed_matrix<double> >
      vcl_ilut(Agpu, ilut_conf);

    xgpu = viennacl::linalg::solve(Agpu,
				   bgpu,
				   custom_bicgstab,
				   vcl_ilut);
  }
  else
    xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab);
  gpuvector2vector(xgpu,x);
  viennacl::backend::finish();
  return x;
}


vector<double> vcl_gmres(sparse::matrix<double>& A, vector<double>& b)
{
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n,0);
  viennacl::vector<double>     bgpu(n), xgpu(n);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::gmres_tag  custom_gmres(1e-7,10000000,50);

  if ( getop("-pre") == "jacobi") {
    viennacl::linalg::jacobi_precond< gpumatrix >
      vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
    xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_gmres, vcl_jacobi);
  }
  else
    xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_gmres);
  gpuvector2vector(xgpu,x);
  viennacl::backend::finish();
  return x;
}
