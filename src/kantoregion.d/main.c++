#include <cstdio>
#include <iostream>
#include <cstring>
#include "est/xmesh.hpp"
using namespace std;


int main()
{
  vector<xyc> Z;
  ifstream ifs;
  ifs.open("kanto.xyc");
  
  in2xyc(ifs,Z);
  ifs.close();

  double miny = 10000000000.0;
  for (int i =1; i< Z.size(); i++) {
    if ( miny> Z[i].y) miny = Z[i].y;
  }
  fprintf(stderr,"miny=%f\n",miny);

  for (int i =1; i<Z.size(); i++) {
    if ( miny == Z[i].y ) Z[i].label = strdup("");
  }

  double maxy = -10000000000.0;
  for (int i =1; i< Z.size(); i++) {
    if ( maxy< Z[i].y) maxy = Z[i].y;
  }
  fprintf(stderr,"maxy=%f\n",maxy);

  for (int i =1; i<Z.size(); i++) {
    if ( maxy == Z[i].y ) Z[i].label = strdup("e2");
  }


  double minx = 10000000000.0;
  for (int i =1; i< Z.size(); i++) {
    if ( minx> Z[i].x) minx = Z[i].x;
  }
  fprintf(stderr,"minx=%f\n",minx);

  for (int i =1; i<Z.size(); i++) {
    if ( abs(minx - Z[i].x)<0.01 ) Z[i].label = strdup("e3");
  }


  double maxx = -10000000000.0;
  for (int i =1; i< Z.size(); i++) {
    if ( maxx < Z[i].x) maxx = Z[i].x;
  }
  fprintf(stderr,"maxx=%f\n",maxx);

  for (int i =1; i<Z.size(); i++) {
    if ( maxx == Z[i].x ) Z[i].label = strdup("");
  }



  for(int i = 1; i<Z.size(); i++){
    if(Z[i].label==NULL) printf("%f %f\n",Z[i].x,Z[i].y);
    else printf("%f %f %s\n",Z[i].x,Z[i].y,Z[i].label);
  }

  double x, y;
  y = miny-0.2;
  for ( x = minx; x< maxx+0.2; x+=0.2)
    if(x==minx) printf("%f %f v0\n",x,y);
    else  printf("%f %f e0\n",x,y);

  
  x = maxx+0.2;
  for ( y = miny-0.2; y<=maxy; y+=0.2)printf("%f %f e1\n",x,y);

  printf("0.16 0.16\n");
  printf("12.3 6.9\n");
  printf("%f %f v3\n",maxx+0.2,maxy);
  printf("12.3 0.16\n");
  return 0;
}
