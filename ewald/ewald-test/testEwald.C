#include <iostream>
#include <vector>
#include "Cell.h"
#include "D3vector.h"
#include "Ewald.h"
#include "Timer.h"
using namespace std;


int main(void)
{
  std::vector<double> vec(9,0);

  vec[0]=4;vec[4]=4;vec[8]=4;

  Cell c(vec);

  Ewald ewald(c);
  Timer tm;
/*
  tm.start();
  ewald.longrange();
  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;
  ewald.print(cout);
*/
  tm.reset();
  tm.start();
  ewald.longrange1();
  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  tm.reset();
  tm.start();
  ewald.shortrange();
  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  tm.reset();
  tm.start();
  ewald.self();
  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;



  ewald.print(cout);


}
