#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "P2_P1.hpp"
#include "mij.hpp"
#include "axij.hpp"
#include "ayij.hpp"
#include "dij.hpp"
#include "hxij.hpp"
#include "hyij.hpp"
#include "est/sparse.hpp"
#include "est/xmesh.hpp"

extern long dimp2(vector<nde>&N);

using namespace sparse;


int main(){
  vector<xyc>Z; vector<nde>N;
  
  f2mesh(fopen("kanto.mesh","r"),Z,N);

  printf("%d\n",dimp2(N));

  
  mij(1,1);



  double u[7],v[7];
  axij(1,1,u,1.0,1.0,1.0);
  ayij(1,1,v,1.0,1.0,1.0);
  dij(1,1);
  hxij(1,1);
  hyij(1,1);
  return 0;
}
