#include <vector>

using namespace std;





class Preconditioner{
public:

  mutable matrix A, L;
  mutable vector<double> d;
  mutable int solver = 0;
  
  
  vector<double>& nsolve(vector<double>& p) const{
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

  vector<double>& icsolve(vector<double>& b) const{
    int n = b.size();
    static vector<double> x(n), y(n);
    L.sync();

    for(int i = 0; i < n; ++i){
      double rly = b[i];
      for ( auto it : L[i] ) { int j = it.first;
	if ( j<i ) rly -= L[i][j]*y[j];
      }
      y[i] = rly/L[i][i];
    }

    for(int i = n-1; i >= 0; --i){
      double lu = 0.0;
      for ( auto it : L[i] ) { int j = it.first;
	if ( i<j && j<n ) lu += L[i][j]*x[j];
      }
      x[i] = y[i]-lu/L[i][i];
    }
    return x;
  }


  vector<double> & solve(vector<double>& b) const {
    switch ( solver ) {
    case 1: return icsolve(b);
    }
    return nsolve(b);
  }

  vector<double>& trans_solve(vector<double>& p) const{
    return solve(p);
  }
  
};
