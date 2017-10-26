#define VIENNACL_WITH_OPENCL
#define CL_DEVICE_DOUBLE_FP_CONFIG
#include <vector>
#include <viennacl/ocl/device.hpp>
#include <viennacl/ocl/utils.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <viennacl/linalg/ichol.hpp>
#include <viennacl/linalg/detail/ilu/chow_patel_ilu.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include "est/sparse.hpp"
#include "cginf.hpp"
//#include "bicgstabinf.hpp"
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
  int n = A.size();
  vector<double> x(n);
  viennacl::compressed_matrix<double> Agpu(n-1,0);
  viennacl::vector<double>     bgpu(n-1), xgpu(n-1);


  printf("CL_DEVICE_DOUBLE_FP_CONFIG=%d\n",CL_DEVICE_DOUBLE_FP_CONFIG);

  if(viennacl :: ocl :: current_device().double_support())
    printf("double support()\n");
  
  matrix2gpumatrix(A,Agpu);
  vector2gpuvector(b,bgpu);

  viennacl::linalg::bicgstab_tag  custom_bicgstab(1e-16,10000,10000);

  custom_bicgstab.abs_tolerance(1e-16);
#if 0
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
#endif
  viennacl::linalg::jacobi_precond< gpumatrix >  vcl_jacobi(Agpu,viennacl::linalg::jacobi_tag());
  
  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab,vcl_jacobi);
  
  viennacl::linalg::bicgstab_tag iters;
  fprintf(stderr,"iter=%d\n",custom_bicgstab.iters());
  fprintf(stderr,"max_iterations=%d\n",custom_bicgstab.max_iterations());
  fprintf(stderr,"max_iterations_before_restart=%d\n",custom_bicgstab.max_iterations_before_restart());
  fprintf(stderr,"error=%f\n",custom_bicgstab.error());
  
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

