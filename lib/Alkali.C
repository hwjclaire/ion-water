#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Alkali.h"
#include "Mlwf.h"
#include "Atom.h"
#include "Cell.h"
#include "D3vector.h"


void Alkali::reset_wf(void)
{ 
  nwf_=0; 
  for ( int i = 0; i < 3; i++) wf_[i] = 0; 
}

void Alkali::reset_atom(void)
{
  xx_ = 0;
  number_ = 0;
}

void Alkali::add_wf( Mlwf & a, bool)
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
  for ( int i = nwf_; i < nwf_; i ++ )
    assert ( wf_[i] == 0 );
#endif

}




void Alkali::add_alkali( Atom & alkali)
{
  //assert ( alkali.mol() == 0 );
  //assert ( xx_ == 0 );
  xx_ = & alkali;
  xx_ -> setmol(this);
}

void Alkali::check_atom ()
{
#if ASSERT
  assert ( xx_ !=0 ); 
  assert ( xx_ -> mol() == this);
#endif
}

void Alkali::check_wf()
{
#if ASSERT
  //cout << nwf_ << endl;
  assert ( nwf_ == 3 );
  for ( int i = 0; i < nwf_; i++ ) 
  {
    assert ( wf_[i] != 0 );
    assert ( wf_[i] -> mol() == this );
  }
#endif
}

double Alkali::distance( Atom* a)
{
  return xx_ -> distance(a,c_);
}

double Alkali::distance( Mlwf* a)
{
  return xx_ -> distance(a,c_);
}

void Alkali::Compute_cm()
{
  cm_ = xx_ -> x();
}  

void Alkali::Compute_cc()
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
    vec = wf_[i] -> x() - xx_ -> x();
    c_->images(vec);
    cc_ -= vec * 2.0;
    charge -= 2.0;
  }

  cc_ /= charge;

  cc_ += xx_ -> x();

}



void Alkali::Compute_dipole()
{
  D3vector vec;
  //cout << cm_ ;
  dipole_ = D3vector(0,0,0);

  tot_charge_ = xx_ -> charge();

  for ( int i = 0; i < nwf_; i ++ )
  {
    vec = wf_[i] -> x() - xx_ -> x();
    c_->images(vec);
    dipole_ -= 2.0 * vec;

  }

  cout << " Warning: Undefined dipole for alkalis\n";

  //cout << nh_ << " " << tot_charge_ << " " << dipole_ << endl;


}  

void Alkali::Compute_polariz()
{

  polariz_.clear();

  for ( int i = 0; i < nwf_; i ++ )
    polariz_ += wf_[i] -> polariz();

}

