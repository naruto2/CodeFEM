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


void squaremesh(int n, vector<xyc> &Z)
{
  double x = 0.0;
  double y = 0.0;
  double h = 0.5;
  
  h = 1/(double)n;
  
  xyc z;

  z.x=0.0,     z.y=0.0,      z.label=(char*)nullptr; Z.push_back(z);
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

      z.x=x, z.y=y, z.label=(char*)nullptr; Z.push_back(z);
   }

  z.x=h/2.0,     z.y=h/2.0,      z.label=(char*)nullptr; Z.push_back(z);
  z.x=1.0-h/2.0, z.y= 1.0-h/2.0, z.label=(char*)nullptr; Z.push_back(z);


}


int main(){
  vector<xyc> Z;
  vector<nde> N;
  squaremesh(6,Z);
  for (unsigned long i=0; i<Z.size(); i++) printf("%ld %f %f %s\n",i,Z[i].x,Z[i].y, Z[i].label);
  int s = Z.size();
  Z.resize(s+3);
  delaunay(Z,N);
  sortmesh(Z,N);
  plotmesh(Z,N);
  return 0;
}
