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
