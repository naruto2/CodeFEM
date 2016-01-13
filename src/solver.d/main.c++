#include <vector>

using namespace std;


class Preconditioner{
public:
  vector<double>& solve(vector<double>& p) const{
    int nn = p.size();
    static vector<double> q(nn);
    for ( int i=0; i<nn; i++ )
      q[i] = p[i];
    return q;
  }

};

#include "solver.hpp"


extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}

int main(int argc, char **argv){
  vector<double> AA(10);
  double B;
  int    JA[10], i, j, n, w;

  genmat(-1,&JA[0],&AA[0],&B);
  n = JA[0]; w = JA[2];

  matrix A(n); vector<double> x(n), b(n);

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
	  A[i-1][JA[j]-1] = AA[j];
	  b[i-1] = B;
	}
      }
  }

  Preconditioner M;
  int max_iter = 100000;
  double double_tol = 0.000001;
  Real tol;

  tol = double_tol;

  for ( i = 0; i< (int)A.size(); i++ ) A[i][0]; //sync
  
  CG(A, x, b, M, max_iter, tol);
  

  
  chkval(stdout,n,&x[0]);
  return 0;
}
