#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Atom.h"
#include "Sulfate.h"
#include "D3vector.h"


void Sulfate::reset_atom(void)
{ 
  s_= h1_= h2_= o_[0] = o_[1] = o_[2] = o_[3] = 0;
  minus1_ = false; minus2_ = false;
}

void Sulfate::reset_wf(void)
{
  nwf_ = 0;
  for ( int i = 0; i < 16; i++) wf_[i] = 0;
}


void Sulfate::add_wf( Mlwf & a)
{
  add_wf( a, true );
}

void Sulfate::add_wf( Mlwf & a, bool noorder)
{
  //cout << nwf_ << endl;
#if ASSERT
  assert ( nwf_ < 16 );
  for ( int i = 0; i < nwf_; i ++) assert( wf_[i] != 0 );
#endif
  wf_[nwf_++] = &a; a.setmol(this);
#if ASSERT
  for ( int i = ++ nwf_; i < 16; i ++) assert( wf_[i] == 0 );
#endif
}

void Sulfate::add_sulfur( Atom & sulfur)
{
  //assert ( sulfur.mol() == 0 );
  //assert ( s_ == 0 );
  s_ = & sulfur;
  s_ -> setmol(this);
}

void Sulfate::add_oxygen( Atom & oxygen )
{
  int no = 0;
  for ( int i = 0; i < 4; i++)
  {
    if ( o_[i] == 0 )
    {
      no = i;
      break;
    }
  }

  o_[no] = & oxygen;
  oxygen.setmol(this);

  //for ( int i = no + 1; i < 4; i++) assert ( o_[i] == 0);
  
      
}
void Sulfate::add_hydrogen( Atom & hydrogen )
{
#if ASSERT
  assert ( hydrogen.mol() == 0 ) ;
  assert ( h2_ == 0 );
#endif
  if ( h1_ == 0 )
  { 
    h1_ = & hydrogen;
  }
  else if( h2_ == 0)
  {
    h2_ = & hydrogen;
  }

  hydrogen.setmol(this);
}



void Sulfate::check_atom ()
{
#if ASSERT
  assert ( s_ !=0 ); 
  assert ( s_ -> mol() == this);
  for ( int i = 0; i < 4; i++)
  {
    assert ( o_[i] != 0 );
    assert ( o_[i] -> mol() == this );
  }

  if ( h1_ == 0 )
  {
    assert( h2_ == 0);
    //minus2_ = true;
  }
  else if ( h2_ == 0)
  {
    //minus1_ = true;
  }
#endif
}

void Sulfate::check_wf()
{
#if ASSERT
  assert ( nwf_ == 16 );
  for ( int i = 0; i < 16; i++ ) 
  {
    assert ( wf_[i] != 0 );
    assert ( wf_[i] -> mol() == this );
  }
#endif
}

double Sulfate::distance( Atom* a )
{
  if ( s_ -> distance ( a, c_ ) > 10) return 10;
  double min = 10;
  for ( int i = 0; i < 4; i++)
  {
    double d = o_[i] -> distance(a,c_);
    if ( d < min ) min = d;
  }
  return min;
}

double Sulfate::distance( Mlwf* a )
{
  if ( s_ -> distance (a,c_) > 10) return 10;
  double min=10;
  for ( int i = 0; i < 4; i++)
  {
    double d = o_[i] -> distance(a,c_);
    if ( d < min ) min = d;
  }
  return min;
}

double Sulfate::distance( Water* a )
{
  return s_ -> distance ( a -> o(), c_);
}

void Sulfate::Compute_cc()
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

  
  for ( int i = 0; i < nwf_; i ++ )
  {
    D3vector vec = wf_[i] -> x() - s_ -> x();

    c_ -> images(vec);

    cc_ -=  vec * 2.0;

    charge -= 2.0;
  }

  cc_ /= charge;

  cc_ += s_ -> x();

}


void Sulfate::Compute_polariz()
{

  polariz_.clear();
  
  for ( int i = 0; i < nwf_; i ++ )
    polariz_ += wf_[i] -> polariz();
  
}

void Sulfate::Compute_cm()//center of mass
{
  cm_[0] = 0; cm_[1] = 0; cm_[2] = 0;

  double mass = s_ -> mass();

  for ( int i = 0; i < 4; i ++ )
  {
    D3vector vec = o_[i] -> x() - s_ -> x();
    
    c_ -> images(vec);  
  
    cm_ += vec * o_[i] -> mass();

    mass += o_[i] -> mass();

    //cout << length(vec) << "  " ;
  }

    //cout << endl;
  
  if ( h2_ != 0 )
  {
    D3vector vec = h2_ -> x() - s_ -> x();

    c_ -> images(vec);

    cm_ += vec * h2_ -> mass();

    mass += h2_ -> mass();
  }

  if ( h1_ != 0 )
  {
    D3vector vec = h1_ -> x() - s_ -> x();

    c_ -> images(vec);

    cm_ += vec * h1_ -> mass();

    mass += h1_ -> mass();
  }
  
  cm_ /= mass;

  cm_ += s_ -> x();

}
  

void Sulfate::Compute_dipole()
{
  D3vector vec = ( s_ -> x() - cm_ ) ;

  tot_charge_ = s_ -> charge();

  c_ -> images(vec);

  dipole_ = vec * s_ -> charge();
    
  for ( int i = 0; i < 4; i ++ )
  {
    D3vector vec = o_[i] -> x() - cm_;

    c_ -> images(vec);  

    dipole_ += o_[i] -> charge() * vec;

    tot_charge_ += o_[i] -> charge();
  }
  


  for ( int i = 0; i < nwf_; i ++ )
  {
    D3vector vec = wf_[i] -> x() - cm_;

    c_ -> images(vec);

    dipole_ -= 2.0 * vec;

    tot_charge_ -= 2.0;
  }

  if ( h2_ != 0 )
  {
    D3vector vec = h2_ -> x() - cm_;

    c_ -> images(vec);

    dipole_ += vec * h2_ -> charge();

    tot_charge_ += h2_ -> charge();
  }

  if ( h1_ != 0 )
  {
    D3vector vec = h1_ -> x() - cm_;

    c_ -> images(vec);

    dipole_ += vec * h1_ -> charge();

    tot_charge_ += h1_ -> charge();
  }

}

void Sulfate::Compute_quad()
{
  quad_.clear();
  quad_ion_.clear();
  quad_mlwfc_.clear();
  quad_mlwfs_.clear();

  int three = 3, ione = 1;
  double one = 1.0;
  char uplo = 'u';
  double charge = s_ -> charge();

  this -> Compute_cc();

  // ion quadrupoles
  D3vector vec = s_ -> x() - cc_;

  c_ -> images(vec);

  dsyr (&uplo, &three, &charge, &vec[0], &ione, &quad_ion_[0], &three);  
  
  for ( int i = 0; i < 4; i ++ )
  {
    D3vector vec = o_[i] -> x() - cc_;

    c_ -> images(vec);

    charge = o_[i] -> charge();

    dsyr (&uplo, &three, &charge, &vec[0], &ione, &quad_ion_[0], &three);
  }

  if ( h2_ != 0 )
  {
    D3vector vec = h2_ -> x() - cc_;

    c_ -> images(vec);

    charge = h2_ -> charge();

    dsyr (&uplo, &three, &charge, &vec[0], &ione, &quad_ion_[0], &three);
  }

  if ( h1_ != 0 )
  {
    D3vector vec = h1_ -> x() - cc_;

    c_ -> images(vec);

    charge = h1_ -> charge();

    dsyr (&uplo, &three, &charge, &vec[0], &ione, &quad_ion_[0], &three);
  }


  // mlwf quadrupole
  for ( int i = 0; i < nwf_; i ++ )
  {
    D3vector vec = wf_[i] -> x() - cc_;

    c_ -> images(vec);
  
    charge = -2;

    dsyr (&uplo, &three, &charge, &vec[0], &ione, &quad_mlwfc_[0], &three);

    // mlwfs quadrupoles
    quad_mlwfs_ += wf_[i] -> quad();
  }
  
  quad_ = quad_ion_ + quad_mlwfc_ + quad_mlwfs_;

  quad_t_ = quad_;
  quad_t_.traceless();

  quad_t_.syev('u', quad_eigval_, quad_eigvec_ ); 
  
}

