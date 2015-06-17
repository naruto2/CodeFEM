#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <stack>
using namespace std;
#include "xyc_nde.h"

static char buf[9999], *S0, S1[999], S2[999], S3[999];
static int forFILE(FILE *fp)
{ 
  if( !feof(fp) ){
    S0 = fgets(buf, 9990, fp);
    if(feof(fp)) return 0;
  }
  S1[0] = '\0';
  S2[0] = '\0';
  S3[0] = '\0';
  sscanf(S0,"%s %s %s",S1,S2,S3);
  return 1;
}
static char *S(int n)
{ 
  if(n==1) return S1[0]=='\0'?NULL:S1;
  if(n==2) return S2[0]=='\0'?NULL:S2;
  if(n==3) return S3[0]=='\0'?NULL:S3;
  if(n==0) return S0;
  return NULL;
}

void fp2xyc(FILE *fp,vector<xyc>&Z)
{ 
  stack<xyc> s;
  xyc x_y_l;
      
  while(forFILE(fp)) if(S(1)!=NULL&&S(2)!=NULL) {
      x_y_l.x     = atof(S(1));
      x_y_l.y     = atof(S(2));
      x_y_l.label =  (S(3)==NULL?NULL:strdup(S(3)));
      s.push(x_y_l);
    }

  Z.resize(s.size()+4);
  
  while(s.size()!=0) {
    Z[s.size()] = s.top();
    s.pop();
  }
}
