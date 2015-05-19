#include <vector>
#include <complex>
#include "D3vector.h"
#include "Cell.h"
#include "Ewald.h"
#include "blas.h"

using namespace std;

Ewald::Ewald ( Cell& cell, int nion, std::vector<D3vector*> pos, 
  std::vector<double*> polar, std::vector<std::vector<double> > & elocal, 
  int nk, double kappa )
: cell_(cell), nkx_(nk), nky_(nk), nkz_(nk), nk_(nkx_*nky_*nkz_-1),
  nion_(nion), kappa_(kappa), ion_(pos), ionpolar_(polar), 
  ionfield_tot_(elocal), pi_ (3.14159265352979), sqrtpi_( 1.77245385091)
{

  

  const double twopi = pi_ * 2;

  const double x = cell_.x();
  const double y = cell_.y();
  const double z = cell_.z();

  const double dkx = twopi / cell_.x() ;
  const double dky = twopi / cell_.y() ;
  const double dkz = twopi / cell_.z() ;
#if 0
  double t;
  cout << "nk?" << endl;
  cin >> t;
  nkx_ = t;
  nky_ = t;
  nkz_ = t;
  nk_ = t * t * t - 1;


  cout << "kappa?" << endl;
  cin >> kappa_;

#endif

  kappa2_ = kappa_ * kappa_;
  kappa3_ = kappa2_ * kappa_;
  fourkappa2_ = 4.0 * kappa2_;
  const_self_ = 4.0 * kappa3_ / 3.0 / sqrtpi_;

  ///////////////////////////////////////////////////////////////////////////
  //Initializa k vectors
  k_.resize(nk_);
  k2_.resize(nk_);
  rho_k_.resize(3*nk_,0);

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
  assert ( ik == nk_);
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  //Initialize local fields
  ionfield_short_.resize(nion_);
  ionfield_long_.resize(nion_);
  ionfield_self_.resize(nion_);

  for ( int i = 0; i < nion_; i ++)
  {
    ionfield_short_[i].resize(9,0);
    ionfield_long_[i].resize(9,0);
    ionfield_self_[i].resize(9,0);
  }
  /////////////////////////////////////////////////////////////////////////////

}


void Ewald::field()
{

  //reset zeros
  for ( int i = 0; i < nion_; i ++)
    for ( int j = 0; j < 9; j ++)
    {
      ionfield_short_[i][j]=0;
      ionfield_long_[i][j]=0;
      ionfield_self_[i][j]=0;
    }

  shortrange();
  longrange();
  self();


  for ( int i = 0; i < nion_; i ++ )
    for ( int j = 0; j < 9; j ++ )
    {
      ionfield_tot_[i][j] = ionfield_short_[i][j] +  ionfield_long_[i][j] + 
                            ionfield_self_[i][j];
    }

  //print(cout);

}




void Ewald::shortrange(void)
{

  int three = 3;
  double one = 1.0;
  char cl = 'l', cn = 'n';
  
  double r2max = cell_.z() * cell_.z();

  #pragma omp parallel
  {
    double eloc[nion_][9];
  
    for ( int iion = 0; iion < nion_ ; iion ++ )
      for ( int i = 0; i < 9; i ++ ) 
        eloc[iion][i]=0;

    #pragma omp for
    for ( int iion = 0; iion < nion_ ; iion ++ ) 
    {
   
      D3vector& r1 = *(ion_ [ iion ]);
   
      for ( int jion = iion + 1; jion < nion_; jion ++ )
      {
   
        D3vector vec0 = *(ion_ [ jion ]) - r1;
    
        cell_.images ( vec0 );
   

        const int ncell = 2;
  
        for ( int ix = - ncell; ix <= ncell; ix ++)
        for ( int iy = - ncell; iy <= ncell; iy ++)
        for ( int iz = - ncell; iz <= ncell; iz ++)
        {

          D3vector vec_cell ( cell_.x() * ix, iy * cell_.y(), iz * cell_.z() );
  
          D3vector vec =  vec0 + vec_cell;
          double r2 = vec * vec;
          if ( r2 > r2max ) continue;
  
          double r = sqrt(r2) , r3 = r2 * r ;
          //double r5 = r2*r3;
     
          //construct T matrix
          double t[9];
          for ( int i = 0; i < 9; i ++ ) t[i] = 0.0;
     
          for ( int i = 0; i < 3; i ++ ) 
            for ( int j = 0; j < 3; j ++ ) 
          {
     
            t[i*3+j] = vec[i] * vec[j] / r2 * ( kappa_ / sqrtpi_ * exp ( -kappa2_ * r2 )
                       * ( 6.0 / r2 + 4.0 * kappa2_ ) + 3.0 / r3 * erfc ( kappa_ * r ) );
     
            if ( i == j )
              t[i*3+j] -= 1.0 / r3 * ( 2 * kappa_ * r / sqrtpi_ * exp ( -kappa2_ * r2 ) 
                          + erfc ( kappa_ * r ) );
     
          }//done T matrix
     
          //cout << iion << "  " << jion << "  " << vec <<  "   " << r << "   " << r2 << endl;
          //for ( int i = 0; i < 9; i ++ )  cout << t[i] << endl;
     
          //add T p_j to E_i and T p_i to E_j
          dgemm ( &cn, &cn, &three, &three, &three, &one, t, &three, ionpolar_[jion], 
            &three, &one,  &eloc[iion][0], &three );
     
          dgemm ( &cn, &cn, &three, &three, &three, &one, t, &three, ionpolar_[iion], 
            &three, &one,  &eloc[jion][0], &three );

        }// for ix, iy, iz
   
      }// for i2
   
    }// for i1

    #pragma omp critical
    for ( int iion = 0; iion < nion_ ; iion ++ )
      for ( int i = 0; i < 9; i ++ )
        ionfield_short_[iion][i]+=eloc[iion][i];


  }

}

void Ewald::longrange(void)
{
  int three = 3, ione = 1;
  char cl = 'l', cn = 'n';
  double one = 1.0, zero = 0.0;
  std::complex<double> imag_unit(0,1);

  const double fourpi = 4 * pi_;
  double fourpioverv = fourpi / cell_.v();

  #pragma omp parallel for
  for ( int ik = 0; ik < 3 * nk_; ik ++ ) 
    rho_k_[ik]=0;

  #pragma omp parallel for
  for ( int ik = 0; ik < nk_; ik ++ )
  {
    D3vector kloc = k_[ik];
    std::complex<double> rhokloc[3]={0,0,0};
    for ( int iion = 0; iion < nion_; iion ++)
    {
      std::complex<double> c1 = imag_unit * exp( kloc * *(ion_[iion]) * imag_unit );

      double* p = &ionpolar_[iion][0];
      for ( int i = 0; i < 3; i ++ )
      {
        complex<double> sum = 0;
        for ( int j = 0; j < 3; j ++ )
          sum += c1 * p[i*3+j] * kloc[j];
        rhokloc[i] += sum;
      }
      //for ( int i = 0; i < 3; i ++ )
      //  rho_k_ [ik*3+i] += c1 * ddot(&three, &ionpolar_[iion][i*3], &ione, &kloc[0], &ione);
    }
    
    #pragma ivdep
    for ( int i = 0; i < 3; i ++ )
    {
      rho_k_ [ik*3+i] = rhokloc[i];
    }
  }

  for ( int iion = 0; iion < nion_; iion ++)
  {

    #pragma omp parallel
    {
      vector <double> fieldloc(9,0);
      

      #pragma omp for
      for ( int ik = 0; ik < nk_; ik ++ )
      {
        double c1 = -fourpioverv * exp ( - k2_[ik] / fourkappa2_ ) / k2_[ik];
        std::complex <double>  ieikr = 
              imag_unit * exp ( k_[ik] * *(ion_[iion]) * imag_unit );
        
        D3vector kloc = k_[ik];

        for ( int i =0; i < 3; i ++ )
        {
          const double fac = c1 * real ( ieikr * conj ( rho_k_[ik*3+i] ) ); 
          double *fieldloci = &fieldloc[i*3];

          //daxpy ( &three, &fac, &k_[ik][0], &ione, &ionfield_long_ [iion][i*3], &ione );
          #pragma ivdep
          for ( int j = 0; j < 3; j ++ )
            fieldloci[j] += kloc[j] * fac;
  
        }
      }

      #pragma omp critical
      for ( int i = 0; i < 9; i ++ )
        ionfield_long_ [iion][i] += fieldloc[i];

    }

  }//for iion

}

void Ewald::self(void)
{
  int nine = 9, ione = 1;
  double zero = 0;

  for ( int iion = 0; iion < nion_; iion ++) 
    daxpy(&nine, &const_self_, ionpolar_[iion], &ione,
            &ionfield_self_ [iion][0], &ione);

}

void Ewald::print(ostream& os)
{
  os.setf(ios::fixed, ios::floatfield);
  os.setf(ios::right, ios::adjustfield);

  os << "Ewald potential Calculation:" << endl
     << " Cell volume: " << cell_.v() << endl
     << " Number of ions: " << nion_ << endl
     << " Number of k vectors: " << nk_ << endl
     << " kappa = " << kappa_ << endl
     << endl;

  for ( int i = 0; i < nion_; i ++ )//mod
  {
    os << i << "   ";
    os << "position: " << setw (15) << *(ion_[i])  << endl;

    os << setw (15) << "polari" 
       << setw (15) << "short" 
       << setw (15) << "long" 
       << setw (15) << "self" 
       << setw (15) << "total" << endl;
            
    for ( int j = 0; j < 9; j ++ )
    {
      os<< setw (15) << ionpolar_[i][j]
        << setw (15) << ionfield_short_[i][j]
        << setw (15) << ionfield_long_[i][j] 
        << setw (15) << ionfield_self_[i][j]
        << setw (15) << ionfield_tot_[i][j] << endl;

    }
  }
}
