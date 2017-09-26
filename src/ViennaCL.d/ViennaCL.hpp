#include <stdio.h>

#ifndef _EST_VIENNACL_HPP_
#define _EST_VIENNACL_HPP_

#define VIENNACL_WITH_OPENCL

#include <iostream>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/row_scaling.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>


using viennacl::linalg::ilut_precond;
using viennacl::linalg::jacobi_precond;
using viennacl::linalg::row_scaling;
using viennacl::compressed_matrix;
using viennacl::linalg::cg_tag;
//using viennacl::linalg::solve;
using viennacl::linalg::ilut_tag;
using viennacl::linalg::jacobi_tag;
using viennacl::linalg::row_scaling_tag;

typedef compressed_matrix<double>  SparseMatrix;
typedef std::vector< std::map< unsigned int, double> > cpumatrix;
//typedef std::vector<double> vector;
typedef viennacl::compressed_matrix<double> gpumatrix;
typedef viennacl::vector<double> gpuvector;
typedef ilut_precond< SparseMatrix >  ILU;
typedef row_scaling< SparseMatrix > Scaling;
typedef jacobi_precond< SparseMatrix > Jacobi;


vector<double> gpucg(matrix& A, vector<double>& b)
{
  int n = A.size();
  vector<double> x(n);
  gpumatrix Agpu(n,0); gpuvector bgpu(n), xgpu(n);
  cpumatrix Acpu(n);

  for (unsigned int i=0; i<A.size(); i++) {
    for ( auto it : A[i] ){
      int j = it.first;
      Acpu[i][j] = A[i][j];
    }
  }
  
  copy(Acpu, Agpu);
  copy(b.begin(), b.end(), bgpu.begin());

  ILU     vcl_ilut(Agpu, ilut_tag(8,1e-3));
  Scaling vcl_row_scaling(Agpu, row_scaling_tag(2));
  Jacobi  vcl_jacobi(Agpu,jacobi_tag());
  cg_tag  custom_cg(1e-5,1000000);

  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_cg, vcl_jacobi);

  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}

#include <viennacl/linalg/bicgstab.hpp>

using viennacl::linalg::bicgstab_tag;

vector<double> gpubicgstab(matrix& A, vector<double>& b)
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
  Jacobi  vcl_jacobi(Agpu,jacobi_tag());
  bicgstab_tag custom_bicgstab(1e-6,100000);

  //viennacl::linalg::ilu0_tag ilu0_config;
  //viennacl::linalg::ilu0_precond< SparseMatrix > vcl_ilut(Agpu, ilu0_config);
  printf("hello\n");
  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_bicgstab, vcl_jacobi);
  printf("world\n");
  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}


#include <viennacl/linalg/gmres.hpp>

using viennacl::linalg::gmres_tag;

vector<double> gpugmres(matrix& A, vector<double>& b)
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
  //Jacobi  vcl_jacobi(Agpu,jacobi_tag());
  //bicgstab_tag custom_bicgstab(1e-12,1000000);
  viennacl::linalg::block_ilu_precond<SparseMatrix, viennacl::linalg::ilu0_tag> vcl_block_ilu0(Agpu, viennacl::linalg::ilu0_tag());
  gmres_tag custom_gmres_tag(1e-8, 100000, 40);
  //viennacl::linalg::ilu0_tag ilu0_config;
  //viennacl::linalg::ilu0_precond< SparseMatrix > vcl_ilut(Agpu, ilu0_config);
  
  xgpu = viennacl::linalg::solve(Agpu, bgpu, custom_gmres_tag, vcl_block_ilu0);
  std::cout << "No. of iters: " << custom_gmres_tag.iters() << std::endl;
  std::cout << "Est. error: " << custom_gmres_tag.error() << std::endl;
  
  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}



#endif


