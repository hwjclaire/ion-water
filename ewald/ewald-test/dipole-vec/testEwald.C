#include <iostream>
#include <vector>
#include <cstdlib>
#include "Cell.h"
#include "D3vector.h"
#include "Ewald.h"
#include "Timer.h"

using namespace std;


int main(int argc, char** argv)
{

  if ( argc != 3 ) return 1;

  double kappa = atof ( argv[2] );
  int nk = atoi ( argv[1] );

  std::vector<double> vec(9,0);

  vec[0]=10;vec[4]=10;vec[8]=10;

  Cell c(vec);

  Ewald ewald(c, nk, kappa);
  Timer tm;

  tm.start();
  ewald.longrange();
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
