////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// D3vector.h
//
// double 3-vectors
//
////////////////////////////////////////////////////////////////////////////////
// $Id: D3vector.h Fri Mar 21 12:09:35 PDT 2014

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "D3vector.h"

using namespace std;

D3vector rotate ( const D3vector& x, const D3vector& w )
{
  if ( length(x) == 0.0 ) return x; // x has zero length
  double theta = length( w );       // rotate by zero
  if ( theta == 0.0 ) return x;
  D3vector ew = normalized ( w );
  D3vector v = w ^ x;
  if ( length( v ) == 0.0 ) return x; // x is parallel to the rotation axis
  v = normalized( v );
  D3vector u = v ^ ew;
  double p = x * u;
  return  (x*ew)*ew + p*cos(theta)*u + p*sin(theta)*v ;
}


double angle ( const D3vector &a, const D3vector&b )
{
  double c = a * b; 
  if ( c > 0 ) 
    return acos(c/length(a)/length(b));
  else
    return - acos(c/length(a)/length(b));
}
  

