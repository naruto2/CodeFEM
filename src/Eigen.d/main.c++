#include "solvers.h"

void fprintmtx(string file, Smatrix &A) {
  FILE *fp;
  if ( file == "stdout" ) fp = stdout;
  else {
    fp = fopen(file.c_str(),"w");
    if( fp == NULL ) return;
  }
  fprintf(fp,"%%%%MatrixMarket matrix coordinate real general\n");

  int n = 0;
  for(int i=0;i<A.outerSize();++i)
    for(Smatrix::InnerIterator it(A,i);it;++it) n++;

  fprintf(fp,"%d %d %d\n", A.rows(), A.cols(), n);

  for(int i=0;i<A.outerSize();++i)
    for(Smatrix::InnerIterator it(A,i);it;++it)
      fprintf(fp,"%d %d %le\n",it.row()+1,it.col()+1,it.value());
  if ( file == "stdout" ) return;
  fclose(fp);
}


int main(void){
#if 1
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
#endif
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

  printmtx("sample0.mtx",A,b);

  Count(ogita);
#endif
  return 0;
}
