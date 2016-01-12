#ifndef _EST_MATRIX_HPP_
#define _EST_MATRIX_HPP_

#include <iostream>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;



class internal_map : public map<int, double> {

  double a;
  int K;

  void wmap(int K, double a) {
    this->erase(K);
    if ( a != 0.0 )
      this->insert( internal_map::value_type(K,a) );
  }

  double rmap( int k) {
    internal_map::iterator it;
    it = this->find(k);
    if ( it == this->end() ) return 0.0;
    return (*it).second;
  }

public:  

  double & operator [](int  k) {
    wmap(K,a);
    K = k;
    a = rmap(k);
    return a;

  }

};


typedef vector< internal_map > internal_matrix;


class matrix : public internal_matrix {

public:

  matrix()      : internal_matrix(){}

  matrix(int n) : internal_matrix(n){}

  friend ostream& operator<<(ostream& os, matrix& A){
    for ( int i=0; i<(int)A.size(); i++ ) {
      for ( int j=0; j<(int)A.size(); j++ ) {
	os << setw(4) << A[i][j] << " ";
      }
      os << endl;
    }
    return os;
  }

#if 0
  friend vector<double>& operator*(matrix A, vector<double>& b) {
    int n = (int)A.size();
    static vector<double> x(n);
    int i, j;
    
    for ( i = 0; i < (int)A.size(); i++) {
      x[i] = 0.0;
      for ( j = 0; j < (int)A.size(); j++) {
	x[i]+=A[i][j]*b[j];
      }
    }
    return x;
  }
#endif  
};

static int dim1(matrix &A){
  return (int)A.size();
}


#endif
