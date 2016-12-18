#include <eigen3/Eigen/Sparse>

#if EIGEN_MAJOR_VERSION<2
#error wrong eigen version
#endif

typedef Eigen::SparseMatrix<double> SpMat;

int main(int argc, char ** argv)
{
  SpMat A(100,100), B;
  B = A.block(10,10,10,10);



}
