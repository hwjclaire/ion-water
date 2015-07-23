#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Water.h"
#include "Mlwf.h"
#include "Atom.h"
#include "Cell.h"
#include "D3vector.h"
#include "blas.h"
//#include "Context.h"
//#include "Matrix.h"

void Water::reset_wf(void)
{ 
  nwf_=0; 
  for ( int i = 0; i < 4; i++) wf_[i] = 0; 
}

void Water::reset_atom(void)
{
  o_ = 0; h_[0] = 0; h_[1] = 0; h_[2] = 0;
  number_ = 0;
}

void Water::add_wf( Mlwf & a, bool noorder)
{
  double thresh = 1;
#if ASSERT
  assert ( nwf_ < 4 );
#endif
  if ( nh_ != 2 || noorder)// wf[0] != 0 && wf[1] != 0 )
  {
    wf_[nwf_] = &a;
    a.setmol(this);
    nwf_++;
    return;
  }
  else
  {

    D3vector oe = a.x() - o_->x();
    c_ -> images(oe);
    //cout << number << length(oe) << endl;
    //oe = normalized(oe);
    //oh1 = normalized(oh1);
    //oh2 = normalized(oh2);
    //cout << oe << "   " << oh1 << "  " <<  oh2;
    //cout  << oh1 * oe <<  "  " << oe * oh2 <<  "  " << oe * ( oh1 + oh2) << endl; 
    //r[nwf] = acos( oh1 * oe ) / 3.1415826 * 180;

    double eoh1 = oe * oh_[0], eoh2 = oe * oh_[1];
    //cout << number_ << "  " << eoh1 << "  " << eoh2 << "  " << length(oe) 
    //     <<  "    " << ( oh1 ^ oh2 ) * oe  <<endl;

    

    if ( eoh1 > thresh ) 
    {
#if ASSERT
      assert ( wf_[0] == 0);
#endif
      wf_[0] = &a;
      a.setmol(this);
      nwf_++;

      return;
    }
    else if ( eoh2 > thresh )
    {
#if ASSERT
      assert ( wf_[1] == 0);
#endif
      wf_[1] = &a;
      a.setmol(this);
      nwf_++;
      return;
    }
    else if ( ( oh_[0] ^ oh_[1] ) * oe > 0 )
    {
#if ASSERT
      assert ( wf_[2] == 0);
#endif
      wf_[2] = &a;
      a.setmol(this);
      nwf_++;
      return;
    }
    else
    {
#if ASSERT
      assert ( wf_[3] == 0);    
#endif
      wf_[3] = &a;
      a.setmol(this);
      nwf_++;
      return;
    }     

  }//end if ( eoh1 && eoh2 ) else
}




void Water::add_oxygen( Atom & oxygen)
{
#if ASSERT  
  assert ( oxygen.mol() == 0 );
  assert ( o_ == 0 );
#endif 
  o_ = & oxygen;
  o_ -> setmol(this);
}

void Water::add_hydrogen( Atom & hydrogen )
{
#if ASSERT  
  assert ( hydrogen.mol() == 0 ) ;
  assert ( nh_ < 3 );
  
  for ( int i = 0; i < nh_; i++)
  {
    assert ( h_[i] != 0);
  }
#endif 
  
  h_[nh_] = & hydrogen;
  hydrogen.setmol(this);

  oh_[nh_] = h_[nh_] ->x() - o_->x();
  c_ -> images(oh_[nh_]);

  nh_ ++;
  
#if ASSERT 
  for ( int i = nh_; i < 3; i++)
  {
    assert ( h_[i] == 0);
  }
#endif

}

void Water::watermove()
{
  for ( int i = 0; i < nh_; i++)
  {
    oh_[i] = h_[i] ->x() - o_->x();
    c_ -> images(oh_[i]);
  }
}



void Water::mlwfmove()
{
  for ( int i = 0; i < nwf_; i++)
  {
    owf_[i] = wf_[i] ->x() - o_->x();
    c_ -> images(owf_[i]);
    
  }

}

void Water::getaxis_dipole()
{
  Compute_dipole();
  
  axis_[0] = normalized ( dipole_ );
  axis_[1] = normalized (oh_[0] ^ oh_[1]);
  axis_[2] = axis_[0] ^ axis_[1];

  for ( int i = 0; i < 3; i ++)
  for ( int j = 0; j < 3; j ++)
    axis_tensor_[i*3+j]=axis_[i][j];


}
void Water::getaxis()
{
  axis_[0] = normalized (oh_[0] + oh_[1]);
  if ( norm(oh_[0]) > norm(oh_[1]) )
    axis_[1] = normalized (oh_[0] ^ oh_[1]);
  else
    axis_[1] = normalized (oh_[1] ^ oh_[0]);
  axis_[2] = axis_[0] ^ axis_[1];
/*
  cout << number_ << endl;
  cout << axis_[0] << endl;
  cout << axis_[1] << endl;
  cout << axis_[2] << endl;
*/
  for ( int i = 0; i < 3; i ++)
  for ( int j = 0; j < 3; j ++)
    axis_tensor_[i*3+j]=axis_[i][j];

  for ( int i = 0; i < nwf_; i ++)
  {

    D3vector v[3];
    v[0] = normalized (owf_[i]);//o-wf vector
    if ( i < 2 )
    {
      v[1] = axis_[1];//vecter perpendicular to the plane
      v[2] = normalized ( v[0] ^ v[1] );
      v[0] = normalized ( v[1] ^ v[2] );

      if ( i % 2 ) v[2] = - v[2];
    }
    else
    {
      v[2] = axis_[2]; // vector in the plane
      v[1] = normalized ( v[0] ^ v[2] );
      v[0] = normalized ( v[2] ^ v[1] );
 
      if ( i % 2 ) v[1] = - v[1];
    }


    for ( int j = 0; j < 3; j ++)
      for ( int k = 0; k < 3; k ++)
        wf_[i]->axis_tensor()[j*3+k]=v[j][k];

    //cout << v[0] << endl;

  }
}


void Water::check_atom ()
{
  /*
  //swap hydrogen so that the first is closet to oxygen.
  if ( o_->distance(h_[0],c_) >  o_->distance(h_[1],c_) )
  {
    Atom * t = h_[1];
    h_[1] = h_[0];
    h_[0] = h_[1];
  }
  */

#if ASSERT  
  assert ( o_ !=0 ); 
  assert ( o_ -> mol() == this);
  assert ( h_[0] !=0 );
  assert ( h_[0] -> mol() == this);
  if( h_[1] == 0 )
  {
    assert( h_[2] == 0);
  }
  else if ( h_[2] == 0 )
  {
    assert ( h_[1] -> mol() == this);
  }
  else
  {
    assert ( h_[2] -> mol() == this);
  }
#endif
}

void Water::check_wf()
{
  //cout << nwf_ << endl;

#if ASSERT
  assert ( nwf_ == 4 );
  for ( int i = 0; i < 4; i++ ) 
  {
    assert ( wf_[i] != 0 );
    assert ( wf_[i] -> mol() == this );
  }
#endif
}

double Water::distance( Atom* a)
{
  return o_ -> distance(a,c_);
}

double Water::distance( Mlwf* a)
{
  return o_ -> distance(a,c_);
}

void Water::Compute_cm()
{
  double mass = 0.0;
  cm_ = D3vector();
  mass = o_->mass();
  for ( int i = 0; i < nh_; i ++ )
  {
    cm_ += oh_[i] * h_[i] -> mass();
    mass += h_[i] -> mass();
  }
  cm_ /= mass;
  
  cm_ += o_ -> x();
}  

void Water::Compute_cc()
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
    cc_ -=  owf_[i] * 2.0;
    charge -= 2.0;
  }

  cc_ /= charge;

  cc_ += o_ -> x();

}

void Water::Compute_polariz()
{

  polariz_.clear();

  for ( int i = 0; i < nwf_; i ++ )
    polariz_ += wf_[i] -> polariz();

}

void Water::Compute_dipole()
{
  //cout << cm_ ;
  dipole_ = D3vector(0,0,0);

  tot_charge_ = o_ -> charge();

  for ( int i = 0; i < nh_; i ++ )
  {
    dipole_ += oh_[i] * h_[i] -> charge();

    //cout << oh_[i] << endl;

    tot_charge_ += h_[i] -> charge();

  }


  for ( int i = 0; i < 4; i ++ )
  {
    dipole_ -= 2.0 * owf_[i];

    //cout << owf_[i] << endl;

    tot_charge_ -= 2.0;
  }

  //cout << nh_ << " " << tot_charge_ << " " << dipole_ <<  " " << length(dipole_) << endl;
  //cout << o_ -> x() << " " << h_[0] -> x() << " " << h_[1] -> x() << endl;

}  

void Water::Compute_quad()
{ 

  quad_.clear();
  quad_ion_.clear();
  quad_mlwfc_.clear();
  quad_mlwfs_.clear();

  int three = 3, ione = 1;
  double one = 1.0;
  double charge = o_ -> charge();
  char uplo = 'u';
  
  //
  {
    D3vector vec = o_->x() - cm_;
    dger (&three, &three, &charge, &vec[0], &ione, &vec[0], &ione, &quad_ion_[0], &three);
  }

  for ( int i = 0; i < nh_; i ++ )
  {
    charge = h_[i] -> charge();

    D3vector vec = oh_[i] + o_->x() - cm_;

    dger (&three, &three, &charge, &vec[0], &ione, &vec[0], &ione, &quad_ion_[0], &three);
  }
 
  for ( int i = 0; i < 4; i ++ )
  {
    //cout << i << endl;
    //cout << owf_[i] << endl;

    charge = -2.0;

    D3vector vec = owf_[i] + o_->x() - cm_;

    dger (&three, &three, &charge, &vec[0], &ione, &vec[0], &ione, &quad_mlwfc_[0], &three);

    quad_mlwfs_ += wf_[i] -> quad() * charge;
    
    //cout << quad_mlwfc_;
    //cout << wf_[i]->quad();

  }

  quad_ = quad_ion_ + quad_mlwfc_ + quad_mlwfs_;

  quad_t_ = quad_;
  quad_t_.traceless();

  quad_t_.syev('u', quad_eigval_, quad_eigvec_ );

}  

void Water::rotate_polar_tensor()
{
  char cn = 'n', ct = 't';
  int three = 3, ione = 1;
  double one = 1.0, zero = 0.0;
  double temp[9];

  dgemm( &cn, &cn, &three, &three, &three, &one, polar_tensor_, &three,
          axis_tensor_, &three, &zero, temp, &three);
  dgemm( &ct, &cn, &three, &three, &three, &one, axis_tensor_, &three,
          temp, &three, &zero, polar_tensor_axis_, &three);
}

void Water::rotate_alpha_tensor()
{
  char cn = 'n', ct = 't';
  int three = 3, ione = 1, nine = 9;
  double one = 1.0, zero = 0.0;
  double temp[9];


#if 1
  dgemm( &cn, &cn, &three, &three, &three, &one, alpha_tensor_, &three,
          axis_tensor_, &three, &zero, temp, &three);
  dgemm( &ct, &cn, &three, &three, &three, &one, axis_tensor_, &three,
          temp, &three, &zero, alpha_tensor_axis_, &three);

#else
  
  Context ctxt_;

  DoubleMatrix a_up(ctxt_,3,3), a_lo(ctxt_,3,3), vec_up(ctxt_,3,3), vec_lo(ctxt_,3,3);
  DoubleMatrix a_sy(ctxt_,3,3), vec_sy(ctxt_,3,3);
  dcopy(&nine, alpha_tensor_, &ione, a_up.valptr(),&ione);
  dcopy(&nine, alpha_tensor_, &ione, a_lo.valptr(),&ione);
  dcopy(&nine, alpha_tensor_, &ione, a_sy.valptr(),&ione);
  std::valarray<double> eig_up, eig_lo, eig_sy;
  eig_up.resize(3,0);
  eig_lo.resize(3,0);
  eig_sy.resize(3,0);

  a_up.syev('u',eig_up,vec_up);
  a_lo.syev('l',eig_lo,vec_lo);
  a_sy.symmetrize('n');
  a_sy.syev('l',eig_sy,vec_sy);

  double test[9];
  dgemm( &ct, &cn, &three, &three, &three, &one, axis_tensor_,
    &three, &vec_sy[0], &three, &zero, test, &three);

  for ( int ii = 0; ii < 3; ii ++ )
  for ( int jj = 0; jj < 3; jj ++ )
  {
    cout <<number_ << "  "<<  ii << "   " << jj << "   " << test [ ii * 3 + jj ];;
    if ( fabs( test [ ii * 3 + jj ] ) > 0.5 )
    {
      cout << "xxxxxxx";
      alpha_tensor_axis_[ii * 4] = eig_sy[jj];
    }
    cout << endl;
  }
#endif

#if 0
#if 0

  for ( int ii = 0; ii < 9; ii ++ )
    cout
       << setw(12) << axis_tensor_[ii]
       << setw(12) << alpha_tensor_[ii]
       << setw(12) << alpha_tensor_axis_[ii]
       << endl;


#else
  for ( int ii = 0; ii < 3; ii ++ ) cout << eig_up[ii] << "   ";
  cout << endl;    
  for ( int ii = 0; ii < 3; ii ++ ) cout << eig_lo[ii] << "   ";
  cout << endl;
  for ( int ii = 0; ii < 3; ii ++ ) cout << eig_sy[ii] << "   ";
  cout << endl;

  for ( int ii = 0; ii < 9; ii ++ )
    cout
       << setw(12) << axis_tensor_[ii]
       << setw(12) << alpha_tensor_[ii]
       << setw(12) << alpha_tensor_axis_[ii]
       << setw(12) << vec_up[ii]
       << setw(12) << vec_lo[ii]
       << setw(12) << vec_sy[ii]
       << setw(12) << test[ii]
       << endl;

#endif
#endif
}

