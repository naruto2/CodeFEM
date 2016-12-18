#include <cstdio>
#include <cmath>
#include <vector>
#include "est/xmesh.hpp"

using namespace std;

static void arrow(FILE *fp,double x0,double y0,double x, double y)
{
  double xl,yl,xr,yr,h;
  h = 0.2;

  xl = (-sqrt(3.0)/2.0)*h; yl = ( 1.0/2.0)*h;
  xr = (-sqrt(3.0)/2.0)*h; yr = (-1.0/2.0)*h;

  fprintf(fp,"%f %f\n",x0,y0);
  fprintf(fp,"%f %f\n",x0+x,y0+y);
  fprintf(fp,"\n");
  fprintf(fp,"%f %f\n",x0+x,y0+y);
  fprintf(fp,"%f %f\n",x0+x+x*xl-y*yl,y0+y+y*xl+x*yl);
  fprintf(fp,"\n");
  fprintf(fp,"%f %f\n",x0+x,y0+y);
  fprintf(fp,"%f %f\n",x0+x+x*xr-y*yr,y0+y+y*xr+x*yr);
  fprintf(fp,"\n");

}



static void plotncpolynomial1_internal(FILE *fp, vector<xyc> Mid, vector<double> u, vector<double> v)
{
  double t, x0,y0,x,y;
  long i,m;
  t = 0.5;

  m = Mid.size()-1;
  for(i=1;i<=m;i++){
    x0 = Mid[i].x;
    y0 = Mid[i].y;

    x  = t*u[i];
    y  = t*v[i];

    arrow(fp,x0,y0,x,y);
    fprintf(fp,"%f %f\n",x0,y0);
    fprintf(fp,"%f %f\n",x0+x,y0+y);
    fprintf(fp,"\n");
  }
  fflush(fp);
}

static const char *filename;
void setanimefilename(const char *fname)
{
  filename=fname;
}

void plotncpolynomial1(vector<xyc> Mid, vector<double> x)
{
  static FILE *pp = NULL;

  if ( filename != NULL && pp == NULL ) {
    char command[512];
    snprintf(command,500,"gzip>%s.gnuplot.gz",filename);
    pp = popen(command,"w");
  }

  if ( pp == NULL ) pp = popen("gnuplot","w");

  unsigned long i, m = Mid.size()-1;
  vector<double> u(m+1), v(m+1);
  for(i=1;i<=m;i++){ u[i] = x[i]; v[i] = x[i+m];}

  fprintf(pp,"plot '-' title \"\" with lines\n");
  plotncpolynomial1_internal(pp, Mid, u, v);
  fprintf(pp,"e\n\n");
  fprintf(pp,"pause 1\n\n");
  fflush(pp);
}

