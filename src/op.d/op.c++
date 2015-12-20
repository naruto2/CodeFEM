#include "op.hpp"

static int    argc;
static char **argv;


void initop(int &pargc, char **&pargv)
{
  argc = pargc; argv = pargv;
}


int defop(string str)
{ 
  for (int i = argc-1; 0<i; i--)
    if ( str == argv[i] )
      return 1;
  return 0;
}


string getop(string str)
{
  int i;

  for ( i=argc-2; 0<i; i-- )
    if ( str == argv[i] )
      return argv[i+1];

  if ( (str == "" ) && 1<argc)
    return argv[i+1];

  return "\0";
}
