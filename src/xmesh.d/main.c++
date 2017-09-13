#include <iostream>
#include "est/xmesh.hpp"
#include "est/op.hpp"

using namespace std;

int main(int argc, char **argv)
{ 
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  ifstream ifs;
  
  if ( defop("-f") ) {
    if (getop("-f") == "-" ) {
      in2xyc(cin,Z);
    }
    else {
      ifs.open(getop("-f"));
      if (!ifs) {
	cerr << "Error: 入力ストリームを開けませんでした(xmesh)" << endl;
	return 0;
      }
      in2xyc(ifs,Z); 
      ifs.close();
    }
  } else {
    in2xyc(cin,Z);
  }

  delaunay(Z, N);
  sortmesh(Z, N);
  p2(Z, N);
  
  if ( defop("-o") ) {
    if ( getop("-o") == "-" ) {
      outmesh(cout,Z,N);
    }
    else {
      ofstream ofs;
      ofs.open(getop("-o").c_str(), ios::out);
      outmesh(ofs,Z,N);
    }
  }
  else
    plotmesh(Z,N);

  return 0;
}
