#ifndef _EST_SPARSE_HPP_
#define _EST_SPARSE_HPP_
#include <iostream>
#include <vector>
#include <map>

using namespace std;

namespace sparse {
  template<typename T>
  class matrix : public vector< map<long, T> > {
  };

  template<typename T>
  void printmatrix(matrix<T>&A)
  {
    long i, j, n;
    n = A.size();
    for(i=1;i<n;i++){
      for(j=1;j<n;j++) cout<<A[i][j]<<" ";
      cout << endl;
    }
  };

}
  
#endif
