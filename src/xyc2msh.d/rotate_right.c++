#include <vector>
using namespace std;
#include "est/xyc2msh.hpp"

void rotate_right(vector<nde>&N,int i)
{ int a,b,c,A,B,C;
  a=N[i].a; b=N[i].b; c=N[i].c; A=N[i].A; B=N[i].B; C=N[i].C;
  N[i].a=c; N[i].b=a; N[i].c=b; N[i].A=C; N[i].B=A; N[i].C=B;
}
