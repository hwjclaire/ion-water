#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Halide.h"
#include "Mlwf.h"
#include "Atom.h"
#include "Cell.h"
#include "D3vector.h"


void Halide::reset_wf(void)
{ 
  nwf_=0; 
  for ( int i = 0; i < 4; i++) wf_[i] = 0; 
}

void Halide::reset_atom(void)
{
  xx_ = 0;
  number_ = 0;
}

void Halide::add_wf( Mlwf & a, bool)
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




void Halide::add_halide( Atom & halide)
{
  //assert ( halide.mol() == 0 );
  //assert ( xx_ == 0 );
  xx_ = & halide;
  xx_ -> setmol(this);
}

void Halide::check_atom ()
{
#if ASSERT
  assert ( xx_ !=0 ); 
  assert ( xx_ -> mol() == this);
#endif
}

void Halide::check_wf()
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

double Halide::distance( Atom* a)
{
  return xx_ -> distance(a,c_);
}

double Halide::distance( Mlwf* a)
{
  return xx_ -> distance(a,c_);
}

void Halide::Compute_cm()
{
  cm_ = xx_ -> x();
}  

void Halide::Compute_cc()
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



void Halide::Compute_dipole()
{
  D3vector vec;
  //cout << cm_ ;
  dipole_ = D3vector(0,0,0);

  tot_charge_ = xx_ -> charge();

  for ( int i = 0; i < 4; i ++ )
  {
    vec = wf_[i] -> x() - xx_ -> x();
    c_->images(vec);
    dipole_ -= 2.0 * vec;

  }

  cout << " Warning: Undefined dipole for halides\n";

  //cout << nh_ << " " << tot_charge_ << " " << dipole_ << endl;


}  

void Halide::Compute_polariz()
{

  polariz_.clear();

  for ( int i = 0; i < nwf_; i ++ )
    polariz_ += wf_[i] -> polariz();

}

