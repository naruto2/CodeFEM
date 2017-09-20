#include <cstdio>
#include "est/xmesh.hpp"

std::vector<xyc> cavity(int n);
  

int main(){
  vector<xyc>Z;
  Z=cavity(8);

  for(int i=1;i<Z.size();i++){
    if (Z[i].label==NULL)
      printf("%f %f\n",Z[i].x,Z[i].y);
    else 
      printf("%f %f %s\n",Z[i].x,Z[i].y,Z[i].label);
  }
}
