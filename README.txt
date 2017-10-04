CodeFEM - Codes for FEM

### 1st.
$ cd bin/; ./mklink; cd ..

### 2nd.
$ cd src/xmesh.d; make

### xmesh is a generator of FEM mesh.

### 3rd.
# mkdir /usr/include/est
# cp include/est/sparse.hpp   /usr/include/est
# cp include/est/bicgstab.hpp /usr/include/est
# cp src/dataparallel.d/cl_bicgstab_kernel.cl /usr/include/est
# cp lib/libbicgstab.a /usr/lib
# c++ sample.c++ -lbicgstab -lOpenCL
# Usage: ./a.out [OPTION]
# argv[1] = Number of the PE. s.t. ./a.out 16

### cl_bicgstab() is a linier solver using OpenCL.
/* sample.c++ --- A sample source for cl_bicgstab(). */
#include <cstdio>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>

int main(int argc, char **argv)
{
  cl_bicgstab_init(argc,argv);

  sparse::matrix<double> A(3);
  vector<double> b(3), x;

  A[1][1] = 2.0; A[1][2] = 1.0; b[1] = 3.0;
                 A[2][2] = 1.0; b[2] = 1.0;

  x = cl_bicgstab(A,b);
  printf("x[1]=%f\n",x[1]);
  printf("x[2]=%f\n",x[2]);
  return 0;
}
