#include <iostream>
#include <iomanip>
#include <vector>
#include "D3vector.h"
#include "Cell.h"
#include "Timer.h"
using namespace std;

int main (void)
{
  std::vector < double > a(9,0), b = a;
  double r = 2.0;
  a[0] = r; a[4] = r; a[8] = r;
  double volume;
  
  Cell c(a);
  
  for ( int i = 0 ; i < a.size(); i ++ )
  {
    //cout << a[i] << endl;
    //cout << b[i] << endl;
  }
  //cout << volume << endl;

  Timer tm;
  tm.start();
  
  int n = 1000000000;
  for ( long int i = 0; i < n; i ++ )  
  {
    D3vector rr(3.212312,312.1239,-1200123.123213123);
    c.images(rr);
    if ( i % (n/10)  == 0) cout << i << endl;
  }

  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;
}
