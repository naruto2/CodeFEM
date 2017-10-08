#include <cstdio>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/navierstokes.hpp>
#include <mcheck.h>

int main(int argc, char **argv)
{
  double Re, dt;

  cl_bicgstab_init(argc,argv);
  if(0!=navierstokes_init("cavity32.mesh",Re=5000,dt=0.001))
    return 0;

  sparse::matrix<double> A; vector<double> U, b;

  for ( int T=0; T<= 36000000; T++) {
    fprintf(stderr,"T");
    navierstokes(A,U,b);

    fprintf(stderr,"=");
    mtrace();
    U = cl_bicgstab(A,b);
    muntrace();
    fprintf(stderr,"%05d\n",T);
    plotuv(U);
    if ( 0 == (T%1000)) fprintuv(U);
  }
  return 0;
}
