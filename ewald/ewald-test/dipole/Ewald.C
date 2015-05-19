#include <vector>
#include <complex>
#include "D3vector.h"
#include "Cell.h"
#include "Ewald.h"
#include "blas.h"

using namespace std;

Ewald::Ewald(Cell& cell, int nk, double kappa) 
  : cell_(cell), nkx_(nk), nky_(nk), nkz_(nk), nk_(nkx_*nky_*nkz_-1)
{
  const double twopi = 2 * 3.14159265358979;

  double x = cell_.x();
  double y = cell_.y();
  double z = cell_.z();

  double dkx = twopi / cell_.x() ;
  double dky = twopi / cell_.y() ;
  double dkz = twopi / cell_.z() ;

  kappa_ = kappa;
  kappa2_ = kappa_ * kappa_;
  kappa3_ = kappa2_ * kappa_;
  fourkappa2_ = 4.0 * kappa2_;

  ///////////////////////////////////////////////////////////////////////////
  //Initializa k vectors
  k_.resize(nk_);
  k2_.resize(nk_);
  rho_k_.resize(nk_,0);

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

  const int nsub = 1;

  nion_ = 2* nsub * nsub * nsub;  

  ion_.resize(nion_);
  ionpolar_.resize(nion_,D3vector(0,0,0));
  ionfield_short_.resize(nion_,D3vector(0,0,0));
  ionfield_long_.resize(nion_,D3vector(0,0,0));
  ionfield_self_.resize(nion_,D3vector(0,0,0));

  int count = 0;
  for ( int i = 0; i < nsub; i ++)
  for ( int j = 0; j < nsub; j ++)
  for ( int k = 0; k < nsub; k ++)
  {
    ionpolar_[count] = D3vector(1,0,0);
    ion_[count] = D3vector (x/nsub*i-.5,y/nsub*j,z/nsub*k);
    count ++;

    ionpolar_[count] = D3vector(1,0,0);
    ion_[count] = D3vector (x/nsub*i+.5,y/nsub*j,z/nsub*k);
    count ++;


    
  }

  assert( count == nion_ );

}

void Ewald::shortrange(void)
{

  int three = 3, ione = 1;
  double one = 1.0, zero = 0.0;
  char cl = 'l', cn = 'n';

  const double fourpi = 2 * 3.14159265358979;
  const double sqrtpi = sqrt(3.14159265358979);

  for ( int iion = 0; iion < nion_ ; iion ++ ) 
    ionfield_short_ [ iion ] = D3vector(0,0,0);

  for ( int iion = 0; iion < nion_ ; iion ++ ) 
  {

    D3vector& r1 = ion_ [ iion ];

    for ( int jion = iion + 1; jion < nion_; jion ++ )
    {

      D3vector vec = ion_ [ jion ] - r1;
  
      cell_.images ( vec );

      double r2 = vec * vec;
      double r = sqrt(r2) , r3 = r2 * r ;
      double r5 = r2*r3;
      if ( r > cell_.x() / 2 ) continue;

      //construct T matrix
      double t[9];
      for ( int i = 0; i < 9; i ++ ) t[i] = 0.0;

      for ( int i = 0; i < 3; i ++ ) 
        for ( int j = i; j < 3; j ++ ) 
        {

#if 1
          t[i*3+j] = vec[i] * vec[j] / r2 * ( kappa_ / sqrtpi * exp ( -kappa2_ * r2 )
                     * ( 6.0 / r2 + 4.0 * kappa2_ ) + 3.0 / r3 * erfc ( kappa_ * r ) );

          if ( i == j )
            t[i*3+j] -= 1.0 / r3 * ( 2 * kappa_ * r / sqrtpi * exp ( -kappa2_ * r2 ) 
                        + erfc ( kappa_ * r ) );

//          if ( jion == iion + 1 )
//          cout << jion << " " << vec[i] << "  " << vec[j] << "   " << vec << "    " << t[i*3+j] << endl;
                
//          t[i*3+j] = vec[i] * vec[j] / r2 * ( kappa_ / sqrtpi * exp ( -kappa2_ * r2 )
//                     * ( 4.0 * kappa2_ ) ) ;
        
#else
          t[i*3+j] = vec[i] * vec[j] / r5 * 3.0;

          if ( i == j )
            t[i*3+j] -= 1.0 / r3;
#endif

        }//done T matrix

      //cout << iion << "  " << jion << "  " << vec <<  "   " << r << "   " << r2 << endl;
      //for ( int i = 0; i < 9; i ++ )  cout << t[i] << endl;

      //add T p_j to E_i and T p_i to E_j
      dsymv ( &cl, &three, &one, t, &three, &ionpolar_[jion][0], &ione, &one, 
                &ionfield_short_ [iion][0], &ione );

      dsymv ( &cl, &three, &one, t, &three, &ionpolar_[iion][0], &ione, &one, 
                &ionfield_short_ [jion][0], &ione );

//      dgemm ( &cn, &cn, &three, &three, &three, &one, t, &three, &ionpolar_[jion][0], 
//        &three, &one,  &ionfield_short_ [iion][0], &three );
/*
      if ( jion  < 16)
      { 
        cout << r1 << endl;
        cout << ion_ [ jion ] << endl;
        cout << vec << endl;
        for ( int i = 0; i < 9; i ++ ) cout << t[i] << " ";
        cout << endl;
      }
*/
    }// for i2

  }// for i1

}

void Ewald::longrange(void)
{
  int three = 3, ione = 1;
  double one = 1.0, zero = 0.0;
  char cl = 'l';
  std::complex<double> imag_unit(0,1);

  const double fourpi = 4 * 3.14159265358979;
  double fourpioverv = fourpi / cell_.v();
  cout << fourpioverv << endl;

  for (int i1 = 0; i1 < nion_ ; i1 ++ ) 
    ionfield_long_ [ i1 ] = D3vector(0,0,0);

  for ( int ik = 0; ik < nk_; ik ++ )
  {
    rho_k_ [ik] = 0;
    for ( int iion = 0; iion < nion_; iion ++)
      rho_k_ [ik] += ionpolar_[iion] * k_[ik] * imag_unit
                      * exp( k_[ik] * ion_[iion] * imag_unit );
  }

  for ( int iion = 0; iion < nion_; iion ++)
  {
  
    for ( int ik = 0; ik < nk_; ik ++ )
    {
      double fac = -fourpioverv * exp ( - k2_[ik] / fourkappa2_ ) / k2_[ik]
                    * real( imag_unit * exp ( k_[ik] * ion_[iion] * imag_unit )
                    * conj ( rho_k_[ik] ) ); 

      daxpy ( &three, &fac, &k_[ik][0], &ione, &ionfield_long_ [iion][0], &ione );
    }

  }//for iion

}



void Ewald::longrange1(void)
{
  int three = 3, ione = 1;
  double one = 1.0, zero = 0.0;
  char cl = 'l';

  const double fourpi = 4 * 3.14159265358979;
  double fourpioverv = fourpi / cell_.v();
  cout << fourpioverv << endl;


  for (int i1 = 0; i1 < nion_ ; i1 ++ ) ionfield_long_ [ i1 ] = D3vector(0,0,0);

  for ( int iion = 0; iion < nion_; iion ++) 
  {
    D3vector& ri = ion_[iion];
    for ( int jion = iion; jion < nion_; jion ++)
    {
      D3vector vec = ion_[jion] - ri;

      //t matrix
      double t[9];
      for ( int i = 0; i < 9; i ++ ) t[i] = 0.0;

      for ( int ik = 0; ik < nk_; ik ++ )
      {

        double fac = -exp ( - k2_[ik] / fourkappa2_ ) / k2_[ik] * cos ( k_[ik] * vec );

        //dger( &three, &three, &fac, &k_[ik][0], &ione, &k_[ik][0], &ione, t, &three );
        dsyr( &cl, &three, &fac, &k_[ik][0], &ione, t, &three );


      } // for k
       
      dsymv ( &cl, &three, &fourpioverv, t, &three, &ionpolar_[jion][0], &ione, &one,
                &ionfield_long_ [iion][0], &ione );

      dsymv ( &cl, &three, &fourpioverv, t, &three, &ionpolar_[iion][0], &ione, &one,
                &ionfield_long_ [jion][0], &ione );

    } //for jion

  }// for iion

  
}

void Ewald::self(void)
{
  const double sqrtpi = sqrt(3.1415926535898);

  for ( int iion = 0; iion < nion_; iion ++) //mod
  {
    ionfield_self_ [iion] = 4.0 * kappa3_ / 3.0 / sqrtpi * ionpolar_[ iion ];
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
     << endl;

  for ( int i = 0; i < nion_; i ++ )//mod
  {
    os << "position: " << setw (15) << ion_[i]  << endl
       << "polari:   " << setw (15) << ionpolar_[i]  << endl
       << "short:    " << setw (15) << ionfield_short_[i]  << endl
       << "long:     " << setw (15) << ionfield_long_[i]  << endl
       << "self:     " << setw (15) << ionfield_self_[i]  << endl
       << "total:    " << setw (15) << 
            ionfield_long_[i] + ionfield_short_[i] + ionfield_self_[i] << endl;

  }
}








