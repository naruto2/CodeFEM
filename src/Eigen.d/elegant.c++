#include <eigen3/Eigen/Sparse>
#include <iostream>

using namespace Eigen;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{

  Matrix4f m = Matrix4f::Random();
  Matrix3f A = Matrix3f::Constant(0.1);
  Vector4f b = Vector4f::Constant(0.2), c = Vector4f::Constant(0.3);
  cout << m << endl << endl;
  cout << A << endl << endl;
  cout << b << endl << endl;
  cout << c << endl << endl;

  m.block(0, 0, 3, 3) = A;
  m.col(3) = b;
  m.row(3) = c;

  cout << m << endl << endl;

  return 0;
}
