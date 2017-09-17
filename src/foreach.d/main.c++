#include <cstdio>
#include <vector>
#include "est/foreach.hpp"
#include "est/xmesh.hpp"
using namespace std;

int main(){
  int i, j;
  double z;
  vector<xyc>Z; vector<nde>N;
  f2mesh(stdin,Z,N);

  for(i=1;i<N.size();i++) {
    printf("%d ",i);
    foreach(j,&N[i].a,&N[i].b,&N[i].c,&N[i].A,&N[i].B,&N[i].C) {
      printf("(");
      foreach(z,&Z[j].x,&Z[j].y) printf("%f ",z);
      printf(")");
      printf("%d ",j);
    }
    printf("\n");
  }
  return 0;
}
