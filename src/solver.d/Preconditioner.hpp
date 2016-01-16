#include <vector>

using namespace std;


class Preconditioner{
public:

  mutable matrix A, L;
  mutable vector<double> d;
  
  vector<double>& nsolve(vector<double>& p) const{
    int n = p.size();
    static vector<double> q(n);
    for ( int i=0; i<n; i++ ) q[i] = p[i];
    return q;
  }

  vector<double>& trans_solve(vector<double>& p) const{
    int n = p.size();
    static vector<double> q(n);
    for ( int i=0; i<n; i++ ) q[i] = p[i];
    return q;
  }

  void jacobi(matrix& pA) {
    A = pA;
  }

  vector<double>& scsolve(vector<double>& p) const{
    int n = p.size();
    static vector<double> q(n);

    for ( int i=0; i<n; i++ ) q[i] = p[i]/A[i][i];
    
    return q;
  }

  vector<double>& solve(vector<double>& p) const{
    int n = p.size();
    static vector<double> q(n);
    static vector<double> y(n);
    backward(L,y,p,d,q);
    return q;
  }

  void setic(matrix Lp, vector<double> dp) {
    L = Lp;
    d = dp;
  }
  
};
