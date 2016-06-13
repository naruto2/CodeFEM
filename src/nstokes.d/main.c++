#include <cstdio>
#include <cstdlib>
#include <cmath>


namespace TaylorHood
{

  double a(int i) {
    switch(i) {
    case 1: return  0.0;
    case 2: return  0.0;
    case 3: return  1.0;
    case 4: return  0.0;
    case 5: return  0.0;
    case 6: return  0.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double b(int i){
    switch(i) {
    case 1: return -1.0;
    case 2: return  0.0;
    case 3: return -3.0;
    case 4: return  0.0;
    case 5: return  4.0;
    case 6: return  0.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double c(int i){
    switch(i) {
    case 1: return  0.0;
    case 2: return -1.0;
    case 3: return -3.0;
    case 4: return  4.0;
    case 5: return  0.0;
    case 6: return  0.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double d(int i){
    switch(i) {
    case 1: return  0.0;
    case 2: return  0.0;
    case 3: return  4.0;
    case 4: return -4.0;
    case 5: return -4.0;
    case 6: return  4.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double e(int i){
    switch(i) {
    case 1: return  2.0;
    case 2: return  0.0;
    case 3: return  2.0;
    case 4: return  0.0;
    case 5: return -4.0;
    case 6: return  0.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double f(int i){
    switch(i) {
    case 1: return  0.0;
    case 2: return  2.0;
    case 3: return  2.0;
    case 4: return -4.0;
    case 5: return  0.0;
    case 6: return  0.0;
    default: abort();
    }
    abort();
    return NAN;
  }
  
  double ad(int j){
    switch(j) {
    case 1: return  0.0;
    case 2: return  0.0;
    case 3: return  1.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double bd(int j){
    switch(j) {
    case 1: return  1.0;
    case 2: return  0.0;
    case 3: return -1.0;
    default: abort();
    }
    abort();
    return NAN;
  }

  double cd(int j){
    switch(j) {
    case 1: return  0.0;
    case 2: return  1.0;
    case 3: return -1.0;
    default: abort();
    }
    abort();
    return NAN;
  }
    
};

using namespace TaylorHood;

int main( int argc, char** argv)
{

  printf("%f\n",a(1));
  return 0;
}

