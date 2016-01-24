#ifndef _EST_PSC98_HPP_
#define _EST_PSC98_HPP_


extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}


void getprob(matrix& A, vector<double>& b);

void check(vector<double>& x);

#endif
