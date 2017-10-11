#ifndef _EST_SPARSE_HPP_
#define _EST_SPARSE_HPP_
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <unistd.h>

using namespace std;

namespace sparse {

  template<typename T>
  class matrix : public vector< map<long, T> > {

  public:
    matrix(){}
    matrix(long n) : vector< map<long, T> >(n){};
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

  template<typename T>
  void plotmatrix(matrix<T>&A)
  {
    FILE *pp;
    long i, j, n=A.size();

    pp = popen("gnuplot","w");
    fprintf(pp,"unset border\n");
    fprintf(pp,"unset xtics\n");
    fprintf(pp,"set xrange [0:%d]\n",n);
    fprintf(pp,"set yrange [%d:0]\n",n);
    fprintf(pp,"set size square\n");
    
    for(i=1;i<n;i++){
      for( auto it : A[i]){
	j = it.first;
	if ( A[i][j] != 0 ) fprintf(pp,"set label \"%.4f\" at %d, %d;\n",A[i][j],j,i);
      }
    }
    fprintf(pp,"plot '-' with lines title \"\"\n");
    fprintf(pp,"1 1\n");
    fprintf(pp,"%d %d\n",n-1,n-1);
    fprintf(pp,"e\n\n");
    fflush(pp);
    sleep(60*3);
    pclose(pp);
  };

}
  
#endif
