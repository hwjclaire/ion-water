#include <vector>
#include <complex>
#include "D3vector.h"
#include "Cell.h"
#include "Ewald.h"

using namespace std;

Ewald::Ewald(Cell& cell) : cell_(cell), nkx_(15), nky_(15), nkz_(15), nk_(nkx_*nky_*nkz_-1)
{
  const double twopi = 2 * 3.14159265358979;
  cout << twopi << endl;

  double x = cell_.x();
  double y = cell_.y();
  double z = cell_.z();

  double dkx = twopi / cell_.x() ;
  double dky = twopi / cell_.y() ;
  double dkz = twopi / cell_.z() ;

  kappa_ = 8/cell_.x();
  fourkappa2_ = 4.0 * kappa_ * kappa_;

  ///////////////////////////////////////////////////////////////////////////
  //Initializa k vectors
  k_.resize(nk_);
  k2_.resize(nk_);
  sk_.resize(nk_,complex<double>(0,0));

  int ik = 0;
  for ( int ikx = -nkx_/2; ikx <= nkx_/2; ikx++)
    for ( int iky = -nky_/2; iky <= nky_/2; iky++)
      for ( int ikz = -nkz_/2; ikz <= nkz_/2; ikz++)
      {
        //cout << ikx << " " << iky << " " << ikz << endl;
        if ( ikx == 0 && iky == 0 && ikz == 0 ) continue;
        k_[ik] = D3vector(ikx*dkx,iky*dky,ikz*dkz);
        //k_[ik] = D3vector(ikx*dkx+dkx/2,iky*dky+dky/2,ikz*dkz+dkz/2);
        k2_[ik] = k_[ik] * k_[ik] ;
        ik ++;
      }
  cout << ik << endl;
  assert ( ik == nk_);
  //////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  //Initializa ions

  const int nsub = 2;


  nion_ =  8 * nsub * nsub * nsub;  

  ion_.resize(nion_);
  ioncharge_.resize(nion_,0);
  ionpotential_short_.resize(nion_,0);
  ionpotential_long_.resize(nion_,0);
  ionpotential_self_.resize(nion_,0);


  //rock salt initialization
  int count = 0;
  for ( int i = 0; i < nsub; i ++)
  for ( int j = 0; j < nsub; j ++)
  for ( int k = 0; k < nsub; k ++)
  {
    ion_[8*count] = D3vector (x/nsub*i,y/nsub*j,z/nsub*k);
    ion_[1+8*count] = D3vector (x/nsub/2+x/nsub*i,0+y/nsub*j,z/nsub/2+z/nsub*k);
    ion_[2+8*count] = D3vector (0+x/nsub*i,y/nsub/2+y/nsub*j,z/nsub/2+z/nsub*k);
    ion_[3+8*count] = D3vector (x/nsub/2+x/nsub*i,y/nsub/2+y/nsub*j,0+z/nsub*k);

    
    ioncharge_[8*count] = 1;
    ioncharge_[1+8*count] = 1;
    ioncharge_[2+8*count] = 1;
    ioncharge_[3+8*count] = 1;


    ion_[4+8*count] = D3vector (x/nsub/2+x/nsub*i,0+y/nsub*j,0+z/nsub*k);
    ion_[5+8*count] = D3vector (0+x/nsub*i,0+y/nsub*j,z/nsub/2+z/nsub*k);
    ion_[6+8*count] = D3vector (0+x/nsub*i,y/nsub/2+y/nsub*j,0+z/nsub*k);
    ion_[7+8*count] = D3vector (x/nsub/2+x/nsub*i,y/nsub/2+y/nsub*j,z/nsub/2+z/nsub*k);

    ioncharge_[4+8*count] = -1;
    ioncharge_[5+8*count] = -1;
    ioncharge_[6+8*count] = -1;
    ioncharge_[7+8*count] = -1;

    count ++;
  }

  assert( count == nion_ / 8);


}

void Ewald::shortrange(void)
{
  // phi^S(r) = 1 / 4pi * sum_j q_j / | r - r_j | * erfc( kappa *  | r - r_j | ) 
  // get phi ( 0 ) first

  const double fourpi = 2 * 3.14159265358979;


  double max,min;

  for ( int i1 = 0; i1 < 1 ; i1 ++ ) //mod // potential at ion i1;
  {
    ionpotential_short_ [ i1 ] = 0;

    D3vector& r1 = ion_ [ i1 ];

    for ( int i2 = 0; i2 < nion_; i2 ++ ) // sum over other ions;
    {
      if ( i1 == i2 ) continue;

      D3vector vec = r1 - ion_ [ i2 ];
  
      cell_.images ( vec );

      double dist = length ( vec );

      double temp = ioncharge_ [ i2 ] / dist * erfc ( dist * kappa_);

      double ftemp = abs ( temp );
    
      if ( i1 == 0 && i2 == 1) { max = ftemp; min = ftemp;}
      //cout << max <<  "  "   << min << endl;

      if ( ftemp < min ) min = ftemp;
      if ( ftemp > max ) max = ftemp;

      ionpotential_short_ [ i1 ] += temp;
  
    }// for i2

  }// for i1


  //cout << max << "    " << min << endl;

}

void Ewald::longrange(void)
{
  // phi^L(k) = sum_j q_j * exp( -i k r_j ) * exp ( k^2 / 4kappa^2 ) / k^2
  // phi^L(r) = sum_{k!=0} exp ( i k r ) * phi^L(k_)

  const double fourpi = 4 * 3.14159265358979;

  for (int i1 = 0; i1 < nion_ ; i1 ++ ) ionpotential_long_ [ i1 ] = 0;

  // get phi^L(k) first
  for ( int ik = 0; ik < nk_; ik ++ )
  {
    for ( int iion = 0; iion < nion_; iion ++)
    {
      double kr = k_[ik] * ion_[iion];
      sk_[ik] += ioncharge_[iion] * exp ( complex<double>(0,-kr) ); 

    }
    //cout << k_ [ik] << "  " << sk_[ik] << endl;
  }

  //get phi^L(r) 
  for ( int iion = 0; iion < 1; iion ++) //mod
  {
    for ( int ik = 0; ik < nk_; ik ++ )
    {
      double kr = k_[ik] * ion_[iion];
      
      ionpotential_long_ [iion] += real ( sk_[ik] / k2_ [ik] * 
                exp ( complex<double>(0,kr) - k2_[ik] / fourkappa2_ ) ) ; 
    }
    ionpotential_long_ [iion] *= fourpi / cell_.v() ;
  }
  
}

void Ewald::longrange1(void)
{
  const double fourpi = 4 * 3.14159265358979;
  double max,min;
  for (int i1 = 0; i1 < nion_ ; i1 ++ ) ionpotential_long_ [ i1 ] = 0;

  for ( int iion = 0; iion < 1; iion ++) //mod
  {
    D3vector& ri = ion_[iion];
    for ( int jion = 0; jion < nion_; jion ++)
    {
      D3vector& rj = ion_[jion];
  
      for ( int ik = 0; ik < nk_; ik ++ )
      {
        //cout << k_[ik]  << "     " << cos ( (ri-rj)*k_[ik] ) << endl;
        double temp = ioncharge_[jion] / k2_[ik] 
            * cos ( (ri-rj)*k_[ik] ) *exp ( - k2_[ik] / fourkappa2_) ;
      
        ionpotential_long_[iion] += temp;

        double t = fabs ( temp * fourpi / cell_.v() );
        
        if ( iion == 0 && jion == 0 && ik == 0 )
        {
          max = t;
          min = t; 
        }
        else
        {
          if ( t > max ) max = t;
          if ( t < min ) min = t;
        }   

      } // for k

    } //for jion

    ionpotential_long_ [iion] *= fourpi / cell_.v();

  }// for iion

  //cout << max << "     " << min << endl;
  
}

void Ewald::self(void)
{
  const double sqrtpi = sqrt(3.1415926535898);

  for (int i1 = 0; i1 < nion_ ; i1 ++ ) ionpotential_self_ [ i1 ] = 0;
 
  for ( int iion = 0; iion < 1; iion ++) //mod
  {
    ionpotential_self_ [iion] -= 2.0 * kappa_ / sqrtpi;
  }

}

void Ewald::print(ostream& os)
{
  const double twopi = 2.0 * 3.14159265358979;


  //os.setf(ios::fixed, ios::floatfield);
  //os.setf(ios::right, ios::adjustfield);

  os << "Ewald potential Calculation:" << endl
     << " Cell volume: " << cell_.v() << endl
     << " Number of ions: " << nion_ << endl
     << " Number of k vectors: " << nk_ << endl
     << " kappa = " << kappa_ << endl
     << " Error Real = " << erfc ( kappa_ * cell_.x() / 2 + 1 ) / ( cell_.x() / 2 + 1 ) << endl
     << " Error Recip = " << exp ( - pow ( twopi / cell_.x() * nkx_ , 2 ) / fourkappa2_ ) 
                              / pow ( twopi / cell_.x() * nkx_ , 2 ) 
     << endl;

  os << setw (12) << "x"
     << setw (12) << "y"
     << setw (12) << "z"
     << setw (15) << "charge"
     << setw (15) << "short"
     << setw (15) << "long"
     << setw (15) << "self"
     << setw (15) << "total" << endl;

  for ( int i = 0; i < 1; i ++ )//mod
  {
    os << setw (15) << ion_[i] 
       << setprecision(10)
       << setw (15) << ioncharge_[i] 
       << setw (15) << ionpotential_short_[i] 
       << setw (15) << ionpotential_long_[i] 
       << setw (15) << ionpotential_self_[i] 
       << setw (15) << ionpotential_long_[i] + ionpotential_short_[i] + ionpotential_self_[i] << endl;

  }
}








