#include <vector>
#include <stack>
using namespace std;
#include "est/xyc2msh.hpp"
#include "incircle.hpp"
#include "Lawson.hpp"

extern void N_set(vector<nde>&,int ,int ,int ,int ,int ,int ,int );
extern int degeneracy(vector<xyc>&,vector<nde>&,int ,int );

void LawsonSwap(vector<xyc>&Z, vector<nde>&N)
{
  int z, i,n;

  z=Z.size()-4;

  N_set(N,1,z+1,z+2,z+3,0,0,0);

  n=1;
  for (i=1;i<=z; i++) {

    int e0,e1,e2,a,b,c,A,B,C; 

    e0=Lawson(Z,N,n,Z[i].x,Z[i].y); e1=n+1; e2=n+2; n+=2;

    a=N[e0].a,b=N[e0].b,c=N[e0].c,A=N[e0].A,B=N[e0].B,C=N[e0].C;

    N_set(N, 0,0,0,0, 0, 0, 0);
    N_set(N,e0,i,b,c,A ,e1,e2);  
    N_set(N,e1,a,i,c,e0,B, e2);
    N_set(N,e2,a,b,i,e0,e1, C);
    
    if(N[B].A==e0) N[B].A=e1;
    if(N[B].B==e0) N[B].B=e1;
    if(N[B].C==e0) N[B].C=e1;

    if(N[C].A==e0) N[C].A=e2;
    if(N[C].B==e0) N[C].B=e2;
    if(N[C].C==e0) N[C].C=e2;
    
    stack<int> s;
    s.push(e0); s.push(e1); s.push(e2);
    
    while(s.size()!=0){
      e1 = s.top();
      s.pop();

      if(incircle(Z,N,i,(e2= (N[e1].a==i?N[e1].A:
			     (N[e1].b==i?N[e1].B:
			      N[e1].C))))){
	if(!degeneracy(Z,N,e1,e2)){ s.push(e1); s.push(e2);}
      }
    }
  }
}
