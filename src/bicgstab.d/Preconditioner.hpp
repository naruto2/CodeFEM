#include <vector>

using namespace std;

class Preconditioner{
public:
  vector<double>& solve(vector<double>& p) const{
    static int i, n = p.size();
    static vector<double> q(n);
    for ( i=1; i<n; i++ ) q[i] = p[i];
    return q;
  }
};
