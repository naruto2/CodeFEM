#ifndef _EST_MATRIX_HPP_
#define _EST_MATRIX_HPP_

#include <iostream>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;



class internal_map : public map<int, double> {

  mutable internal_map *m;
  mutable int K;
  mutable double a;
  mutable internal_map::iterator it;
  
  void wmap(int Key, double a) const {
    m = const_cast<internal_map*>(this);
    m->erase(Key);
    if ( a != 0.0 ) m->insert( internal_map::value_type(Key,a) );
  }

  double rmap(int k) const {
    m = const_cast<internal_map*>(this);
    it = m->find(k);
    if ( it == m->end() ) return 0.0;
    return it->second;
  }

public:  

  double & operator[](const int  k) const {
    wmap(K,a);
    K = k;
    a = rmap(k);
    return a;
  }

};


typedef vector< internal_map > internal_matrix;


class matrix : public internal_matrix {

public:

  internal_matrix AT;
  
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

  void sync() {
    for ( int i =0; i<(int)this->size(); i++ ) (*this)[i][0];
  }
  
  void T() {
    int n = this->size();
    AT.resize(n);
    internal_map aaa;
    internal_map::iterator it;

    for ( int i=0;i<n;i++){
      aaa = (*this)[i];
      it = aaa.begin();
      while(it != aaa.end()){
	AT[it->first][i] = it->second;
	it++;
      }
    }
    for ( int i =0; i<n; i++ ) AT[i][0];
  }

  vector<double>& trans_mult( vector<double> &y) const {
    int n = y.size();
    static vector<double> x(n);

    for ( int i=0; i<n; i++ ) {
      internal_map ATi = AT[i];
      internal_map::iterator j = ATi.begin();
      for ( x[i]=0.0; j != ATi.end(); j++ )
	x[i] += j->second * y[j->first];
    }
    return x;
  }
  
};

static int dim1(matrix &A){
  return (int)A.size();
}


#endif
