#include <vector>

using namespace std;


class Preconditioner{
public:

  matrix A;

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

  vector<double>& solve(vector<double>& p) const{
    int n = p.size();
    static vector<double> q(n);

    for ( int i=0; i<n; i++ ) q[i] = p[i]/A[i][i];
    
    return q;
  }
  
};
