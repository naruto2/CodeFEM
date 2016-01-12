#include "est/xmesh.hpp"
#include "est/op.hpp"


int main(int argc, char **argv)
{ 
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  ifstream ifs(argv[1]);

  if (!ifs) {
    cerr << "Error: 入力ストリームを開けませんでした(xmesh)" << endl;
    return 0;
  }

  in2xyc(ifs,Z); 
  ifs.close();

  delaunay(Z, N);

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
