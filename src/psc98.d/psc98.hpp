#ifndef _EST_PSC98_HPP_
#define _EST_PSC98_HPP_

extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}

//void getprob(matrix& A, vector<double>& b);
void check(vector<double>& x);

void psc98_init(sparse::matrix<double>&A,vector<double>&b);
void psc98_check(vector<double>&x);
#endif
