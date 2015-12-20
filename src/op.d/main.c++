#include <iostream>

#include "est/op.hpp"


int main(int argc, char **argv){

  initop(argc,argv);

  if (defop("-v"))
    cout << getop("-v") << " " << getop("-f") << endl;

  return 0;
}
