#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <stack>
#include <iostream>
#include <fstream>
using namespace std;
#include "xyc_nde.h"

static char S1[999], S2[999], S3[999];

static int forFILE(ifstream &ifs)
{ 
  string str;
  
  if (!getline(ifs,str)) return 0;
  if ( str.length() > 990 ) {
    fprintf(stderr,
	    "Error: In fp2xyc(). The length of the line is too long.\n");
    abort();
  }
  S1[0] = '\0'; S2[0] = '\0'; S3[0] = '\0';
  sscanf(str.c_str(),"%s %s %s",S1,S2,S3);
  return 1;
}

static char *S(int n)
{ 
  if(n==1) return S1[0]=='\0'?NULL:S1;
  if(n==2) return S2[0]=='\0'?NULL:S2;
  if(n==3) return S3[0]=='\0'?NULL:S3;
  return NULL;
}

void ifs2xyc(ifstream &ifs,vector<xyc>&Z)
{ 
  stack<xyc> s;
  xyc x_y_l;

  while(forFILE(ifs)) if(S(1)!=NULL&&S(2)!=NULL) {
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
