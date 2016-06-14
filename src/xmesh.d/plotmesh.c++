#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <unistd.h>
#include "est/xmesh.hpp"

#define  dim1(N) (long)(N.size()-1)

using namespace std;


void plotmesh(vector<xyc>&Z,vector<nde>&N)
{
  FILE *fp = stdout;
  unsigned long e, v, a, b, c;

  fp = popen("gnuplot","w");

  for(e=1;e<N.size();e++){
    fprintf(fp,"set label \"%ld\" at %f, %f;\n",e,
	    (Z[N[e].a].x+Z[N[e].b].x+Z[N[e].c].x)/3.0,
	    (Z[N[e].a].y+Z[N[e].b].y+Z[N[e].c].y)/3.0);
  }

  for(v=1;v<Z.size();v++){
    fprintf(fp,"set label \"%ld\" at %f, %f;\n",v,
	    Z[v].x,
	    Z[v].y);
  }

  fprintf(fp,"plot '-' title \"\" with lines\n");
  
  for(e=1;e<N.size();e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"%f %f\n",Z[b].x,Z[b].y);
    fprintf(fp,"%f %f\n",Z[c].x,Z[c].y);
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"\n\n");
  }

  fprintf(fp,"e\n");
  fflush(fp);
  sleep(65535);
  pclose(fp);
  exit(0);
}
