#include "est/sparse.hpp"

extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}


static void getprob(sparse::matrix<double>& A, vector<double>& b) {
  vector<double> AA(10);
  double B;
  int    JA[10], i, j, n, w;

  genmat(-1,&JA[0],&AA[0],&B);
  n = JA[0]; w = JA[2];

  A.resize(n+1);  b.resize(n+1);

  for ( i=1; i<=n; i++) {
    for (j=0;j<=w-1;j++) {
      JA[j] =  -1;
      AA[j] = 0.0;
    }
    genmat(i,&JA[0],&AA[0],&B);
    for ( j=0; j<w; j++) if (JA[j] != -1){
	if( JA[j] <=0 || n < JA[j] ) {
	  ;
	} else {
	  A[i][JA[j]] = AA[j];
	  b[i] = B;
	}
      }
  }
}

void psc98_init(sparse::matrix<double>& A, vector<double>& b)
{
  getprob(A,b);
}

static void check(vector<double>& x) {
  chkval(stdout,x.size()-1,&x[1]);
}

void psc98_check(vector<double>& x)
{
  chkval(stdout,x.size()-1,&x[1]);
}
