
#include <iostream>

using namespace std;
#include "Context.h"
#include "Matrix.h"

int main ()
{
  Context ctxt;

  DoubleMatrix a(ctxt,3,3);
  DoubleMatrix b(ctxt,3,3);
  DoubleMatrix c(ctxt,3,3);
  DoubleMatrix d(b);


  a.identity();
  b.identity();
  
  a.set('u',0.5);
  a.set('l',0.5);
  b.set('l',1);

  c.gemm('n','n',1.0,a,b,0.0);

  a.inverse();


  d.gemm('n','n',1.0,c,b,0.0);


  cout << a << b << c << d;

  
  double aa[9];
  for ( int i = 0; i < 9; i ++ ) aa[i]=i;

  d.init( aa, 3 );

  cout << d;

}
