#include <stdio.h>

#ifndef _EST_VIENNACL_HPP_
#define _EST_VIENNACL_HPP_

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
using viennacl::linalg::solve;
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


vector<double> gpusolver(matrix& A, vector<double>& b)
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
  cg_tag  custom_cg(1e-7,1000000);

  xgpu = solve(Agpu, bgpu, custom_cg, vcl_jacobi);

  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}

#endif
