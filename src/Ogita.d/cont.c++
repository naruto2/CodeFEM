#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "est/op.hpp"

using namespace std;

void contop(int &argc0, char** &argv0){
  string in_str;
  ifstream ifs(".cont");
  getline(ifs, in_str);

  stringstream ss(in_str);
  char delim = ' ';
  string str;
  static vector<string> v;
  while(getline(ss, str, delim)){
    v.push_back(str);
  };

  char **argv = (char**)malloc(sizeof(char*)*v.size());
  int argc = v.size();
  for ( unsigned int i = 0; i < v.size(); i++){
    argv[i] = (char*)v[i].c_str();
  }
  initop(argc,argv);
  argc0 = argc; argv0 = argv;
}


void mkcont(int argc, char**argv)
{
  FILE *fp = fopen(".cont","w");
  for(int i=0; i< argc; i++) {
    fprintf(fp,"%s ",argv[i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
}


#if 0
int main(int argc, char** argv)
{
  initop(argc,argv);
  mkcont(argc,argv);
  contop();
  return 0;
}
#endif
