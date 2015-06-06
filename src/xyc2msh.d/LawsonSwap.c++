#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"

#define T double
extern int incircle(vector<xyc>&, vector<nde>&, int ,int);
extern void N_set(vector<nde>&,int ,int ,int ,int ,int ,int ,int );
extern int Lawson(vector<xyc>&, vector<nde>&,int , T , T );
extern int degeneracy(vector<xyc>&,vector<nde>&,int ,int );

static int Q[100000], top=0;

void pushq(int e)
{
  Q[top++]=e;
  if ( top > 99997 ) { fprintf(stderr,"LowsonSwap Stack over flow\n");exit(0);}
}

void *popq(int *e)
{
  *e = Q[--top];
  if (top < 0 ) {top = 0; return NULL;}
  return e;
}

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
    
    pushq(e0); pushq(e1); pushq(e2); 

    while(popq(&e1)){

      if(incircle(Z,N,i,(e2=(N[e1].a==i?N[e1].A:(N[e1].b==i?N[e1].B:N[e1].C))))){
	if(!degeneracy(Z,N,e1,e2)){ pushq(e1); pushq(e2);}
      }
    }
  }
}
