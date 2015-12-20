#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#include "op.hpp"
#include "xyc_nde.h"


extern void in2xyc(istream &, vector<xyc>&);
extern void delaunay(vector<xyc>&, vector<nde>&);
extern void outmesh(ostream &, vector<xyc>&, vector<nde>&);
extern void plotmesh(vector<xyc>&, vector<nde>&);


int main(int argc, char **argv)
{ 
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  ifstream ifs(argv[1]);

  if (!ifs) {
    cerr << "Error: 入力ストリームを開けませんでした(xyc2msh)" << endl;
    return 0;
  }

  in2xyc(ifs,Z); 
  ifs.close();

  delaunay(Z, N);

  if ( defop("-o") ) {
    ofstream ofs;
    ofs.open(getop("-o").c_str(), ios::app);
    outmesh(ofs,Z,N);
  }
  else
    plotmesh(Z,N);

  return 0;
}
