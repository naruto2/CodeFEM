#ifndef _EST_MATRIX_HPP_
#define _EST_MATRIX_HPP_

#include <cstdio>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;


class internal_map : public map<int, double> {

  double a;
  int K;

public:
  void wmap(int K, double a) {
    this->erase(K);
    if ( a != 0.0 )
      this->insert( internal_map::value_type(K,a) );
  }

  double rmap( int k) {
    internal_map::iterator it;
    it = this->find(k);
    return (*it).second;
  }
  
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
    for ( int i=0; i<A.size(); i++ ) {
      for ( int j=0; j<A.size(); j++ ) {
	os << setw(4) << A[i][j] << " ";
      }
      os << endl;
    }
    return os;
  }

};



#endif
