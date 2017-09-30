#include <stdio.h>
#include <vector>
#include "est/sparse.hpp"

#define forall(m,i,n) for(i=m;i<=n;i++)
#define less(x,y)     (x<y?(x):(y))
#define absv(x)       ((x)>0.0?(x):-(x))


int LU(sparse::matrix<double>& a)
{  double nrm,aik,apkk;

  static long i,j,k,w,n,pk,pklimit,klimit;
  n=a.size()-1; w=n;

  vector<long>p(n+1);

  for(k=1;k<=n;k++){
    pk=k; nrm=absv(a[k][k]);

    forall(k+1,i,less(k+w,n))if(nrm<absv(a[i][k])){
      pk=i; nrm=absv(a[i][k]);
    }
    if(nrm == 0.0) return 2;
    p[k]=pk;
    klimit=less(k+w,n); apkk= -a[pk][k];
    forall(k+1,i,klimit){
      a[i][k] /= apkk;
      pklimit=less(pk+w,n); aik=a[i][k]; auto apk=a[pk]; auto ai=a[i];
      forall(k+1,j,pklimit) ai[j] += aik*apk[j];
    }
  }
  return 0;
}
