#include <vector>
using namespace std;
#include "xyc2msh.hpp"

void rotate_left(vector<nde>&N,int i)
{ int a,b,c,A,B,C;
  a=N[i].a; b=N[i].b; c=N[i].c; A=N[i].A; B=N[i].B; C=N[i].C;
  N[i].a=b; N[i].b=c; N[i].c=a; N[i].A=B; N[i].B=C; N[i].C=A;
}
