#include <iostream>
#include <time.h>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;

#include "EigenTools.h"
#include "bicgstab.h"
#include "cg.h"
#include "ogita.h"
#include "cgs.h"
#include "bicg.h"
#include "qmr.h"
#include "Vsolvers.h"
#include "psc98.h"


int main(void){
  Smatrix A; Vector b;
  psc98(A,b);
  printf("n = %d\n",b.size());

  Count(Vcg);  
  Count(Vbicgstab);
  Count(Vgmres);
  Count(cg);
  Count(cgs);
  Count(bicg);
  Count(bicgstab);
  Count(qmr);
  Count(Elu);
  Count(Ecg);
  Count(Ebicgstab);
  return 0;

#if 0
  const double s=19, u=21, p=16, e=5, r=18, l=12;
  Tri(0,0, s);
  Tri(1,0, l);
  Tri(4,0, l);
  Tri(1,1, u);
  Tri(2,1, l);
  Tri(4,1, l);
  Tri(0,2, u);
  Tri(2,2, p);
  Tri(0,3, u);
  Tri(3,3, e);
  Tri(3,4, u);
  Tri(4,4, r);
  Smatrix A   = MapSmatrix(5,5);
  double  B[] = {1,1,1,1,1};
  Vector  b   = MapVector(B, 5);
  Count(ogita);
#endif
  return 0;
}
