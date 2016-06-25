#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <cmath>
#include "est/op.hpp"
#include "est/xmesh.hpp"

using namespace std;


static int eq(double a, double b)
{
  if ( fabs(a-b) < 0.01) return 1;
  return 0;
}


int main(int argc, char **argv) {
  initop(argc,argv);
  double x = 0.0;
  double y = 0.0;
  double h = 0.5;
  int n;
  
  if(defop("-n")) n = atoi(getop("-n").c_str());

  h = 1/(double)n;
  
#if 0
  switch (n) {
  case 2 : h = 0.5     ; break;
  case 4 : h = 0.25    ; break;
  case 8 : h = 0.125   ; break;
  case 16: h = 0.0625   ; break;
  case 32: h = 0.03125 ; break;    
  case 64: h = 0.015625  ; break;
  default: h = 0.5; break;
  }
#endif
  
  vector<xyc> Z;
  xyc z;

  
  for ( x = 0.0; x < 1.0+h/2.0; x += h )
    for ( y = 0.0; y < 1.0+h/2.0; y += h ) {
      if ( eq(x , 0.0) && eq(y , 0.0) ) { z.x=x, z.y=y, z.label=(char*)"v0"; Z.push_back(z); continue; }
      if ( eq(y , 0.0) && eq(x , 1.0) ) { z.x=x, z.y=y, z.label=(char*)"v1"; Z.push_back(z); continue; }
      if ( eq(y , 1.0) && eq(x , 1.0) ) { z.x=x, z.y=y, z.label=(char*)"v2"; Z.push_back(z); continue; }
      if ( eq(y , 1.0) && eq(x , 0.0) ) { z.x=x, z.y=y, z.label=(char*)"v3"; Z.push_back(z); continue; }

      if ( eq(y , 0.0) ) { z.x=x, z.y=y, z.label=(char*)"e0"; Z.push_back(z); continue; }
      if ( eq(x , 1.0) ) { z.x=x, z.y=y, z.label=(char*)"e1"; Z.push_back(z); continue; }
      if ( eq(y , 1.0) ) { z.x=x, z.y=y, z.label=(char*)"e2"; Z.push_back(z); continue; }
      if ( eq(x , 0.0) ) { z.x=x, z.y=y, z.label=(char*)"e3"; Z.push_back(z); continue; }

      z.x=x, z.y=y, z.label=(char*)""; Z.push_back(z);
   }

  z.x=h/2.0,     z.y=h/2.0,      z.label=(char*)""; Z.push_back(z);
  z.x=1.0-h/2.0, z.y= 1.0-h/2.0, z.label=(char*)""; Z.push_back(z);

  for(unsigned long i = 0; i<Z.size(); i++)
    printf("%f %f %s\n",Z[i].x, Z[i].y, Z[i].label);
  return 0;
}
