#include <future>
#include <vector>
#include <iostream>
#include <unistd.h>

using namespace std;

bool test(int x, future<int> &f) {
  for(unsigned int i=0; i<100000; i++ ){
    printf("hello %d\n",i);
    f.get();
  }
  return true;
}

int main(int argc, char **argv) {

  promise<int> p;
  future<int> f = p.get_future();

  future<bool> fut = async(launch::async, test, 10, ref(f));

  sleep(3);
  p.set_value(8);

  sleep(3);
  cout << fut.get() << endl;

  return 0;
}
