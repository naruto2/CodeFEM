#include <iostream>
#include <time.h>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> Smatrix;
typedef VectorXd Vector;
void Tri(int i, int j, double v);
Smatrix MapSmatrix(int m, int n);
Vector MapVector(double *array, int n);
Smatrix catSmatrix(Smatrix A, Smatrix B);
Smatrix vcatSmatrix(Smatrix C, Smatrix D);
Vector catVector(Vector x, Vector y);
Vector halfVector(Vector xx);
int IsSymmetric(Smatrix &A);

// Preconditioner M; 前処理クラス
class Preconditioner{
 public:

  Vector solve(Vector p) const{
    Vector q;
    q = p;
    return q;
  }

  Vector trans_solve(Vector p) const{
    Vector q;
    q = p;
    return q;
  }
};

double norm(const Vector &b);
double dot(const Vector &a, const Vector &b);
double check(Vector method(Smatrix &A,Vector &b), Smatrix &A, Vector &b);
#define Count(method) count(#method, method, A, b)
void
count( string str, Vector method(Smatrix &A,Vector &b), Smatrix &A, Vector &b);
Vector Elu(Smatrix &A, Vector &b);
Vector Ebicgstab(Smatrix &A, Vector &b);
Vector Ecg(Smatrix &A, Vector &b);
Vector Eogita(Smatrix &A, Vector &b);

Vector bicgstab(Smatrix &A, Vector &b);
Vector cg(Smatrix &A, Vector &b);
Vector ogita(Smatrix &A, Vector &b);
Vector cgs(Smatrix &A, Vector &b);
Vector bicg(Smatrix &A, Vector &b);
Vector qmr(Smatrix &A, Vector &b);

#define VIENNACL_WITH_OPENCL
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/gmres.hpp>
typedef viennacl::compressed_matrix<double> gpumatrix;
typedef viennacl::vector<double> gpuvector;
typedef std::vector< std::map< unsigned int, double> > cpumatrix;
using namespace viennacl::linalg;
gpuvector btogpu(Vector &b);
Vector fromgpu(gpuvector xgpu);
gpumatrix Atogpu(Smatrix &A);
Vector Vbicgstab(Smatrix &A, Vector &b);
Vector Vcg(Smatrix &A, Vector &b);
Vector Vgmres(Smatrix &A, Vector &b);
Vector Vogita(Smatrix &A, Vector &b);

extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}
void psc98(Smatrix &A, Vector &b);
void checkpsc98(Vector &x);

#include "Esolvers.h"
#include "Vsolvers.h"
#include "Dsolvers.h"
#include "psc98.h"
#include "MatrixMarket.h"
