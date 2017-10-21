/* =========================================================================
   Copyright (c) 2010-2016, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.
                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------
   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
   (A list of authors and contributors can be found in the PDF manual)
   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */
typedef double ScalarType;
#include "viennacl/vector.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/tools/matrix_generation.hpp"
#include "viennacl/linalg/amg.hpp"
#include <iostream>
#include <vector>
#include <ctime>
#include "vector-io.hpp"
#include "viennacl/tools/timer.hpp"
using namespace std;

template<typename MatrixType, typename VectorType, typename SolverTag, typename PrecondTag>

viennacl::vector<double> run_solver(MatrixType const & matrix, VectorType const & rhs, VectorType const & ref_result, SolverTag const & solver, PrecondTag const & precond)
{
  VectorType result(rhs);
  VectorType residual(rhs);
  viennacl::tools::timer timer;
  timer.start();
  result = viennacl::linalg::solve(matrix, rhs, solver, precond);
  viennacl::backend::finish();
  std::cout << "  > Solver time: " << timer.get() << std::endl;
  return result;
  residual -= viennacl::linalg::prod(matrix, result);
  std::cout << "  > Relative residual: " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(rhs) << std::endl;
  std::cout << "  > Iterations: " << solver.iters() << std::endl;
  result -= ref_result;
  std::cout << "  > Relative deviation from result: " << viennacl::linalg::norm_2(result) / viennacl::linalg::norm_2(ref_result) << std::endl;
}
template<typename ScalarType>
viennacl::vector<double> run_amg(viennacl::linalg::cg_tag & cg_solver,
             viennacl::vector<ScalarType> & vcl_vec,
             viennacl::vector<ScalarType> & vcl_result,
             viennacl::compressed_matrix<ScalarType> & vcl_compressed_matrix,
             std::string info,
             viennacl::linalg::amg_tag & amg_tag)
{
  std::cout << "-- CG with AMG preconditioner, " << info << " --" << std::endl;
  viennacl::linalg::amg_precond<viennacl::compressed_matrix<ScalarType> > vcl_amg(vcl_compressed_matrix, amg_tag);
  std::cout << " * Setup phase (ViennaCL types)..." << std::endl;
  viennacl::tools::timer timer;
  timer.start();
  vcl_amg.setup();
  viennacl::backend::finish();
  std::cout << "  > Setup time: " << timer.get() << std::endl;
  std::cout << " * CG solver (ViennaCL types)..." << std::endl;
  return run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, vcl_amg);
}
#include <vector>
using namespace std;
#include "est/sparse.hpp"
#include "est/psc98.hpp"
#include "togpu.hpp"

int main(int argc, char **argv)
{
  std::string filename("mat65k.mtx");
  if (argc == 2)
    filename = argv[1];
  std::cout << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "               Device Info" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
#ifdef VIENNACL_WITH_OPENCL
  // Optional: Customize OpenCL backend
  viennacl::ocl::platform pf = viennacl::ocl::get_platforms()[0];
  std::vector<viennacl::ocl::device> const & devices = pf.devices();
  // Optional: Set first device to first context:
  viennacl::ocl::setup_context(0, devices[0]);
  // Optional: Set second device for second context (use the same device for the second context if only one device available):
  if (devices.size() > 1)
    viennacl::ocl::setup_context(1, devices[1]);
  else
    viennacl::ocl::setup_context(1, devices[0]);
  std::cout << viennacl::ocl::current_device().info() << std::endl;
  viennacl::context ctx(viennacl::ocl::get_context(0));
#else
  viennacl::context ctx;
#endif
  typedef double    ScalarType;  // feel free to change this to double if supported by your device
  viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(ctx);
  //viennacl::tools::generate_fdm_laplace(vcl_compressed_matrix, points_per_dim, points_per_dim);
  // Read matrix
  std::cout << "Reading matrix..." << std::endl;
  std::vector< std::map<unsigned int, ScalarType> > read_in_matrix;
  if (!viennacl::io::read_matrix_market_file(read_in_matrix, filename))
  {
    std::cout << "Error reading Matrix file" << std::endl;
    return EXIT_FAILURE;
  }
  sparse::matrix<double> A; vector<double> b, x;
  //viennacl::copy(read_in_matrix, vcl_compressed_matrix);
  psc98_init(A,b); x.resize(A.size());
  A[0][0] = 1.0;
  matrix2gpumatrix(A,vcl_compressed_matrix);
  std::cout << "Reading matrix completed." << std::endl;
  viennacl::vector<ScalarType> vcl_vec(vcl_compressed_matrix.size1(), ctx);
  viennacl::vector<ScalarType> vcl_result(vcl_compressed_matrix.size1(), ctx);
  std::vector<ScalarType> std_vec, std_result;
  // rhs and result vector:
  std_vec.resize(vcl_compressed_matrix.size1());
  std_result.resize(vcl_compressed_matrix.size1());
  for (std::size_t i=0; i<std_result.size(); ++i)
    std_result[i] = ScalarType(1);
  // Copy to GPU

  //viennacl::copy(std_vec, vcl_vec);
  viennacl::copy(std_result, vcl_result);
  //vcl_vec = viennacl::linalg::prod(vcl_compressed_matrix, vcl_result);


  //vector2gpuvector(b,vcl_vec);
  vector2gpuvector(b,vcl_vec);
  
  viennacl::linalg::cg_tag cg_solver(1e-8, 10000);
  viennacl::context host_ctx(viennacl::MAIN_MEMORY);
  viennacl::context target_ctx = viennacl::traits::context(vcl_compressed_matrix);
  std::cout << "-- CG solver (no preconditioner, warmup) --" << std::endl;
  viennacl::vector<double> xgpu;
  xgpu = run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, viennacl::linalg::no_precond());

  gpuvector2vector(xgpu,x);
  psc98_check(x);


  viennacl::linalg::amg_tag amg_tag_direct;
#if 0
  amg_tag_direct.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_ONEPASS);
  amg_tag_direct.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_DIRECT);
  amg_tag_direct.set_strong_connection_threshold(0.25);
  amg_tag_direct.set_jacobi_weight(0.67);
  amg_tag_direct.set_presmooth_steps(1);
  amg_tag_direct.set_postsmooth_steps(1);
  amg_tag_direct.set_setup_context(host_ctx);    // run setup on host
  amg_tag_direct.set_target_context(target_ctx); // run solver cycles on device
  run_amg(cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "ONEPASS COARSENING, DIRECT INTERPOLATION", amg_tag_direct);
#endif


  
  viennacl::linalg::amg_tag amg_tag_agg_pmis;
  amg_tag_agg_pmis.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
  amg_tag_agg_pmis.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_AGGREGATION);
  xgpu = run_amg(cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), AG INTERPOLATION", amg_tag_agg_pmis);

  gpuvector2vector(xgpu,x); psc98_check(x);



  viennacl::linalg::amg_tag amg_tag_sa_pmis;
#if 0  
  amg_tag_sa_pmis.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
  amg_tag_sa_pmis.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_SMOOTHED_AGGREGATION);
  xgpu = run_amg (cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);
#endif

  std::cout << std::endl;
  std::cout << " -------------- Benchmark runs -------------- " << std::endl;
  std::cout << std::endl;
  std::cout << "-- CG solver (no preconditioner) --" << std::endl;
  xgpu = run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, viennacl::linalg::no_precond());
  gpuvector2vector(xgpu,x); psc98_check(x);
  
  xgpu = run_amg(cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "ONEPASS COARSENING, DIRECT INTERPOLATION", amg_tag_direct);
  gpuvector2vector(xgpu,x); psc98_check(x);
  
  xgpu = run_amg(cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), AG INTERPOLATION", amg_tag_agg_pmis);
  gpuvector2vector(xgpu,x); psc98_check(x);

  
  xgpu = run_amg (cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);
  gpuvector2vector(xgpu,x); psc98_check(x);
  
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
  return EXIT_SUCCESS;
}

