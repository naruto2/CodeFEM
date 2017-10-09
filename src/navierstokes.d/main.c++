#include <cstdio>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/navierstokes.hpp>

double u(double x, double y)
{
  if ( islabel("v0")||islabel("v1")||islabel("v2")||islabel("v3") )
    return 0.0;

  if ( islabel("e0")||islabel("e1")||islabel("e3") )
    return 0.0;

  if ( islabel("e2") )
    return 1.0;

  fprintf(stderr,"A wrong label is exist.\n");
  abort();
}


double v(double x, double y)
{
  if ( islabel("v0")||islabel("v1")||islabel("v2")||islabel("v3") )
    return 0.0;

  if ( islabel("e0")||islabel("e1")||islabel("e2")||islabel("e3") )
    return 0.0;
  
  fprintf(stderr,"A wrong label is exist.\n");
  abort();
}


int main(int argc, char **argv)
{
  double Re, dt;

  cl_bicgstab_init(argc,argv);
  if(0!=navierstokes_init("cavity32.mesh",Re=5000,dt=0.001, u, v))
    return 0;

  sparse::matrix<double> A; vector<double> U, b;

  for ( int T=0; T<= 36000000; T++) {
    fprintf(stderr,"T");
    navierstokes(A,U,b);

    fprintf(stderr,"=");
    U = cl_bicgstab(A,b);

    fprintf(stderr,"%05d\n",T);
    plotuv(U);
    if ( 0 == (T%1000)) fprintuv(U);
  }
  return 0;
}
