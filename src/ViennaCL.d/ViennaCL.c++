#define VIENNACL_WITH_OPENCL
#define CL_DEVICE_DOUBLE_FP_CONFIG
#include <vector>
#include <viennacl/ocl/device.hpp>
#include <viennacl/ocl/utils.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <viennacl/linalg/ichol.hpp>
#include <viennacl/linalg/detail/ilu/chow_patel_ilu.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include "est/sparse.hpp"
//#include "cginf.hpp"
//#include "bicgstabinf.hpp"
#include "gmresinf.hpp"
#include "togpu.hpp"
#include "est/op.hpp"
#include "est/TDMA.hpp"

int isSymmetric(sparse::matrix<double>& A)
{
  int n = A.size();
  for ( int i = 1; i<n; i++) 
    for ( auto it: A[i] ) {
      int j = it.first;
      if ( it.second != A[j][i] ) return 0;
    }
  return 1;
}


vector<double> vcl_cg(sparse::matrix<double>& A, vector<double>& b)
{
  int n = A.size(), diag = 1;
  vector<double> x(n);

  for ( int k=1, n=A.size(); k<n; k++)
    if (A[k][k] == 0.0) { diag = 0; break; }
  
  if ( diag && isTridiagonal(A) ) return TDMA(A,b);

  if ( !isSymmetric(A) ) {
    fprintf(stderr,"Warning: vcl_cg() can't solve asymmetic matrix\n");
    return x;
  }

  viennacl::compressed_matrix<double> Agpu(n-1,0);
  viennacl::vector<double>     bgpu(n-1), xgpu(n-1);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::jacobi_precond< gpumatrix >
    vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());

  viennacl::linalg::cg_tag  custom_cg(1e-14,A.size()/8);
  
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
  int n = A.size(), diag = 1;

  for ( int k=1, n=A.size(); k<n; k++)
    if (A[k][k] == 0.0) { diag = 0; break; }
  
  if ( diag && isTridiagonal(A) ) return TDMA(A,b);

  viennacl::compressed_matrix<double> Agpu(n-1,0);
  viennacl::vector<double>   bgpu(n-1), xgpu(n-1);
  
  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::bicgstab_tag  custom_bicgstab(1e-15,A.size()/8);

  if (diag) {
    viennacl::linalg::jacobi_precond< gpumatrix >
      vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
    xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab,vcl_jacobi);
  } else
    xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab);

  vector<double> x(n);
  gpuvector2vector(xgpu,x);

  if ( defop("-v" ) ) {
    //fprintf(stderr,"iter=%d ",custom_bicgstab.iters());
    //fprintf(stderr,"max_iterations=%d\n",custom_bicgstab.max_iterations());
    }
  viennacl::backend::finish();
  return x;
}



vector<double> vcl_gmres(sparse::matrix<double>& A, vector<double>& b)
{
  int n = A.size(), diag=1;
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n-1,0);
  viennacl::vector<double>     bgpu(n-1), xgpu(n-1);

  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::gmres_tag  custom_gmres(1e-8,A.size()/8,50);

  if ( diag ) {
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


using viennacl::linalg::bicgstab_tag;
using viennacl::linalg::jacobi_precond;
typedef jacobi_precond< sparse::matrix<double> > Jacobi;


vector<double> gpubicgstab(sparse::matrix<double>& A, vector<double>& b)
{

  int n = A.size();
  vector<double> x(n);
  gpumatrix Agpu(n,0); gpuvector bgpu(n), xgpu(n);
  cpumatrix Acpu(n);

  Acpu.clear();
  Acpu.resize(n);


  for (unsigned int i=0; i<A.size(); i++) {
    for ( auto it : A[i] ){
      int j = it.first;
      Acpu[i][j] = A[i][j];
    }
  }

  copy(Acpu, Agpu);
  copy(b.begin(), b.end(), bgpu.begin());
  //ILU          vcl_ilut(Agpu, ilut_tag(256,1e-14));
  //Scaling vcl_row_scaling(Agpu, row_scaling_tag(2));
  bicgstab_tag custom_bicgstab(1e-12,1000000);

  //viennacl::linalg::ilu0_tag ilu0_config;
  //viennacl::linalg::ilu0_precond< SparseMatrix > vcl_ilut(Agpu, ilu0_config);

  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab);


  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}

