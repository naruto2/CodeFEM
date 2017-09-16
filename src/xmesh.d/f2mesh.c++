#include <vector>
#include <cstdio>
#include <cstring>
#include "est/xmesh.hpp"

using namespace std;

static FILE *fcopy(void){
  int c;
  FILE *fp = tmpfile();

  while ( EOF != ( c = getchar() )) putc(c,fp);
  rewind(fp);
  return fp;
}

void f2mesh(FILE *fp, vector<xyc>&Z, vector<nde>&N){
  int m=0, n=0, i;
  char *buf, *label;
  double x, y;
  int a,b,c,A,B,C;

  if ( fp == stdin ) fp = fcopy();

  buf   = (char*)calloc(sizeof(char),1000);
  label = (char*)calloc(sizeof(char),1000);

  while(strncmp(buf,"<xyc>",5)) { fgets(buf,999,fp); }

  while(1) {
    fgets(buf,999,fp);
    if(!strncmp(buf,"<nde>",5)) break;
    sscanf(buf,"%d",&i);
    n = max(n,i);
  } 
  while(1) {
    fgets(buf,999,fp);
    if (feof(fp)) break;
    sscanf(buf,"%d",&i);
    m = max(m,i);
  } 
  rewind(fp);
  Z.resize(n+1); N.resize(m+1);

  while(strncmp(buf,"<xyc>",5)) { fgets(buf,999,fp); }

  while(1) {
    fgets(buf,999,fp);
    if(!strncmp(buf,"<nde>",5)) break;
    label[0] = 0;
    sscanf(buf,"%d %lf %lf %s",&i,&x,&y,label);
    Z[i].label = strdup(label);
    Z[i].x = x;
    Z[i].y = y;
  } 
  while(1) {
    fgets(buf,999,fp);
    if (feof(fp)) break;
    sscanf(buf,"%d %ld %ld %ld %ld %ld %ld",&i,
	   &a,&b,&c,
	   &A,&B,&C);
    N[i].a = a;
    N[i].b = b;
    N[i].c = c;
    N[i].A = A;
    N[i].B = B;
    N[i].C = C;
  } 
  free(buf);
  free(label);
  fclose(fp);
}
