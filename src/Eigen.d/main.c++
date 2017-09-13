#include "solvers.h"

int main(int argc, char **argv){
  Smatrix A; Vector b;
  fscanmtx(argv[1],A,b);

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
}
