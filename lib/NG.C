#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "NG.h"
#include "Mlwf.h"
#include "Atom.h"
#include "Cell.h"
#include "D3vector.h"


void NG::reset_wf(void)
{ 
  nwf_=0; 
  for ( int i = 0; i < 4; i++) wf_[i] = 0; 
}

void NG::reset_atom(void)
{
  o_ = 0;
  number_ = 0;
}

void NG::add_wf( Mlwf & a)
{

#if ASSERT
  for ( int i = 0; i < nwf_; i ++ )
    assert ( wf_[i] != 0 );
#endif

  wf_[nwf_] = &a;
  a.setmol(this);
  a.setnumber(nwf_);
  nwf_++;

#if ASSERT
  for ( int i = nwf_; i < 4; i ++ )
    assert ( wf_[i] == 0 );
#endif

}




void NG::add_oxygen( Atom & oxygen)
{
  //assert ( oxygen.mol() == 0 );
  //assert ( o_ == 0 );
  o_ = & oxygen;
  o_ -> setmol(this);
}

void NG::check_atom ()
{
#if ASSERT
  assert ( o_ !=0 ); 
  assert ( o_ -> mol() == this);
#endif
}

void NG::check_wf()
{
#if ASSERT
  //cout << nwf_ << endl;
  assert ( nwf_ == 4 );
  for ( int i = 0; i < 4; i++ ) 
  {
    assert ( wf_[i] != 0 );
    assert ( wf_[i] -> mol() == this );
  }
#endif
}

double NG::distance( Atom* a)
{
  return o_ -> distance(a,c_);
}

double NG::distance( Mlwf* a)
{
  return o_ -> distance(a,c_);
}

void NG::Compute_cm()
{
  cm_ = o_ -> x();
}  

void NG::Compute_cc()
{
  double charge = 0.0;
  cc_ = D3vector();

/*  
  for ( int i = 0; i < nh_; i ++ )
  {
    cc_ += oh_[i] * h_[i] -> charge();

    charge += h_[i] -> charge();
  }
*/
  D3vector vec;

  for ( int i = 0; i < nwf_; i ++ )
  {
    vec = wf_[i] -> x() - o_ -> x();
    c_->images(vec);
    cc_ -= vec * 2.0;
    charge -= 2.0;
  }

  cc_ /= charge;

  cc_ += o_ -> x();

}



void NG::Compute_dipole()
{
  D3vector vec;
  //cout << cm_ ;
  dipole_ = D3vector(0,0,0);

  tot_charge_ = o_ -> charge();

  for ( int i = 0; i < 4; i ++ )
  {
    vec = wf_[i] -> x() - o_ -> x();
    c_->images(vec);
    dipole_ -= 2.0 * vec;

  }

  //cout << nh_ << " " << tot_charge_ << " " << dipole_ << endl;


}  
