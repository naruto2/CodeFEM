#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
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
static void FILE_cp(FILE *fp,FILE *out)
{ char c; 
  while(EOF !=(c=getc(fp))) putc(c,out);
}     

void fp2xyc(FILE *fp,vector<xyc>&Zv)
{ int i,z; FILE *tfp;

  FILE_cp(fp,(tfp=tmpfile()));
  rewind(tfp);

  z=0;
  while(forFILE(tfp)) if(S(1)!=NULL&&S(2)!=NULL) z++;
  rewind(tfp);

  Zv.resize(z+4);

  i= 1;
  while(forFILE(tfp)) if(S(1)!=NULL&&S(2)!=NULL) {
      Zv[i].x = atof(S(1)); Zv[i].y = atof(S(2));
      Zv[i].label =  (S(3)==NULL?NULL:strdup(S(3)));
      i++;
    }
  fclose(tfp);
}
