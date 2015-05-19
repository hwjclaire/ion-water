#include <algorithm>
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include "D3vector.h"
#include "Mol.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "blas.h"
#include "Context.h"
#include "Matrix.h"
#include "Polariz.h"
#include "Ewald.h"

using namespace std;


Polariz::Polariz(Cell& c, std::vector < Mol* > & molset, 
  std::vector < Mlwf > & mlwfset)
: c_(c), molset_(molset), mlwfset_(mlwfset),
  nmol_(molset.size()), nmlwf_(mlwfset.size()),
  elocal_inv_(ctxt_,3,3)
{

  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::right, ios::adjustfield);
  cout.precision(8);
 

  // need more control for conervergence
  // this parameters are pretty arbitrary
  // may leads to bad convergence for some systems
  nk_ = 21;
  kappa_ = .12;


  //initialize vectors for ewald
  mol_elocal_.resize(nmol_);
  mol_pos_.resize(nmol_);
  mol_polar_.resize(nmol_);


  mlwf_elocal_.resize(nmlwf_);
  mlwf_pos_.resize(nmlwf_);
  mlwf_polar_.resize(nmlwf_);


  for ( int iw = 0; iw < nmol_; iw ++ )
  {
    mol_elocal_ [iw].resize(9);
    mol_pos_ [iw] = &(molset_[iw]->cc());//important
    mol_polar_ [iw] = molset_[iw]->polar_tensor();
    //cout << iw << " " << * mol_pos_ [iw] << endl;
  }

  for ( int i = 0; i < nmlwf_; i ++ )
  {
    mlwf_elocal_ [i].resize(9);
    mlwf_pos_[i] = &(mlwfset_[i].x());
    mlwf_polar_[i] = mlwfset_[i].polar_tensor();
  }


  ewald_mol_ = new Ewald ( c_, nmol_, mol_pos_, mol_polar_, 
                                mol_elocal_, nk_, kappa_ );

  ewald_mlwf_ = new Ewald ( c_, nmlwf_, mlwf_pos_, mlwf_polar_,
                                mlwf_elocal_, nk_, kappa_ );
  
  //intialize external field
  for ( int i = 0; i < 9; i ++ ) eext_[i] = 0;
  eext_[0] = 1.0;
  eext_[4] = 1.0;
  eext_[8] = 1.0;

}

void Polariz::Calculate_mlwf_Ewald(void)
{
  int nine = 9, ione = 1, three = 3;
  double one = 1.0, zero = 0.0;
  char cn = 'n', ct = 't';

  Timer tm;

  //calculate ewald to get ion contibution to the local field
  //these contributions are stored in mlwf_elocal_ and
  //add with external field to get total local field
  tm.start();
  ewald_mlwf_ -> field();
  tm.stop();
  //cout << "Ewald: Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  tm.reset();
  tm.start();
  //calculate total local field
  for ( int i = 0; i < nmlwf_; i ++ )
  {

    //local field
    double * elocal = &mlwf_elocal_[i][0]; 
    daxpy ( & nine, &one, eext_, &ione, elocal, &ione );

    //inverse local field
    elocal_inv_.init(elocal,3);

    elocal_inv_.inverse();


    //alpha * elocal = polarization
    //alpha = polarization * elocal^{-1}
    dgemm( &cn, &cn, &three, &three, &three, &one, mlwf_polar_[i], &three,
          elocal_inv_.valptr(), &three, &zero, mlwfset_[i].alpha_tensor(), &three);

    //cout << "Calc done\n";

    //Rotation
    double temp1[9];
    dgemm( &cn, &cn, &three, &three, &three, &one, mlwfset_[i].alpha_tensor(), &three,
          mlwfset_[i].axis_tensor(), &three, &zero, temp1, &three);
    dgemm( &ct, &cn, &three, &three, &three, &one, mlwfset_[i].axis_tensor(), &three,
          temp1, &three, &zero, mlwfset_[i].alpha_tensor_axis(), &three);

#if 1

    //test
    double test[9];
    dgemm( &cn, &cn, &three, &three, &three, &one, mlwfset_[i].alpha_tensor(),
      &three, elocal, &three, &zero, test, &three);

    cout << i << "    " << mlwfset_[i].number() << endl;
    for ( int ii = 0; ii < 9; ii ++ )
      cout << setw(14) << eext_[ii]
           << setw(14) << elocal[ii]
           << setw(14) << elocal_inv_[ii]
           << setw(14) << mlwfset_[i].alpha_tensor()[ii]
           << setw(14) << mlwfset_[i].alpha_tensor_axis()[ii]
           << setw(14) << mlwfset_[i].axis_tensor()[ii]
           << setw(14) << mlwf_polar_[i][ii]
           << setw(14) << test[ii] << endl;
#endif

  }//for i
  tm.stop();
  //cout << "Other: Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;
}

void Polariz::Polar_mlwf_rotate(void)
{
  int nine = 9, ione = 1, three = 3;
  double one = 1.0, zero = 0.0;
  char cn = 'n', ct = 't';

  //calculate total local field
  for ( int i = 0; i < nmlwf_; i ++ )
  {

    //Rotation
    double temp1[9];
    dgemm( &cn, &cn, &three, &three, &three, &one, mlwfset_[i].polar_tensor(), &three,
          mlwfset_[i].axis_tensor(), &three, &zero, temp1, &three);
    dgemm( &ct, &cn, &three, &three, &three, &one, mlwfset_[i].axis_tensor(), &three,
          temp1, &three, &zero, mlwfset_[i].polar_tensor_axis(), &three);

#if 0

    cout << i << "    " << mlwfset_[i].number() << endl;
    for ( int ii = 0; ii < 9; ii ++ )
      cout 
           << setw(14) << mlwfset_[i].polar_tensor()[ii]
           << setw(14) << mlwfset_[i].polar_tensor_axis()[ii]
           << setw(14) << mlwfset_[i].axis_tensor()[ii] << endl;
#endif

  }//for i
}

void Polariz::Calculate_mol_Ewald(void)
{
  int nine = 9, ione = 1, three = 3;
  double one = 1.0, zero = 0.0;
  char cn = 'n', ct = 't';

  Timer tm;

  //initialize mol polarizations
  for ( int iw = 0; iw < nmol_; iw ++ )
  { 
    for ( int i = 0; i < 9; i ++)
      mol_polar_[iw][i] = 0;

    for ( int imlwf = 0; imlwf < molset_[iw]->nwf(); imlwf ++ )
      daxpy ( &nine, &one, molset_[iw]->wf(imlwf) -> polar_tensor(),
                &ione, mol_polar_[iw], &ione );
  }

  //calculate ewald to get ion contibution to the local field
  //these contributions are stored in mol_elocal_ and
  //add with external field to get total local field
  tm.start();
  ewald_mol_ -> field();
  tm.stop();
  //cout << "Ewald: Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  tm.reset();
  tm.start();
  //calculate total local field
  for ( int i = 0; i < nmol_; i ++ )
  {

    //local field
    double * elocal = &mol_elocal_[i][0]; 
    daxpy ( & nine, &one, eext_, &ione, elocal, &ione );

    //copy w to w_inverse and inverse
    dcopy(&nine, elocal, &ione, elocal_inv_.valptr(),&ione);
    elocal_inv_.inverse();

    //alpha * elocal = polarization
    //alpha = polarization * elocal^{-1}
    dgemm( &cn, &cn, &three, &three, &three, &one, mol_polar_[i], &three,
          elocal_inv_.valptr(), &three, &zero, molset_[i]->alpha_tensor(), &three);

    //cout << "Calc done\n";

  }//for i
  tm.stop();
  //cout << "Other: Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;
}



//Lattice sum plus 4 pi / 3 * P 

void Polariz::Calculate_mol(void)
{

  int nine = 9, ione = 1, three = 3; 
  double one = 1.0, zero = 0.0;
  char cn = 'n', ct = 't';
  double fourpi = 4.0 * 3.1415926535898 / c_.v();
  //depolarization ration 4 pi / 3
  double depolar = 4.0 * 3.14159 / 3.0 / c_.v();

  //external field = maxwell field + 4 pi polarization
  //E0 = E + 4pi*P
  for ( int imol = 0; imol < nmol_; imol ++ )
  for ( int imlwf = 0; imlwf < molset_[imol]->nwf(); imlwf ++ )
  {
    daxpy ( &nine, &depolar, molset_[imol]->wf(imlwf) -> polar_tensor(),
                &ione, eext_, &ione );
  }

  for ( int i = 0; i < nmol_; i ++ )//mod
  {
    //cout << i << endl;
    D3vector moli = molset_[i]->cc();

    //cout << molset_[i].wf(0)->polar()[0] << endl;
    
    //effective poarizability tensor of mol i
    double polari[9];
    for ( int ii = 0; ii < 9; ii ++ ) polari[ii] = 0.0;

    for ( int imlwf = 0; imlwf < molset_[i]->nwf(); imlwf ++ )
    {
      daxpy ( &nine, &one, molset_[i]->wf(imlwf) -> polar_tensor(),
                &ione, &polari[0], &ione );

    }

    //elocal store local field of WF i
    double  elocal[9];
    for ( int ii = 0; ii < 9; ii ++ ) elocal[ii] = 0.0;
    
    #pragma omp parallel
    {
    double  elocal_loc[9];
    for ( int ii = 0; ii < 9; ii ++ ) elocal_loc[ii] = 0.0;

    #pragma omp for
    for ( int j = 0; j < nmol_; j ++ ) //mod 
    {

      const D3vector& molj = molset_[j]->cc();
      //cout << wfj << endl;


      //effective polarizability tensor of mol j
      double polarj[9];
      for ( int ii = 0; ii < 9; ii ++ ) polarj[ii] = 0.0;
      

      for ( int imlwf = 0; imlwf < molset_[j]->nwf(); imlwf ++ )
        daxpy ( & nine, &one, molset_[j]->wf(imlwf) -> polar_tensor(), 
                &ione, polarj, &ione );
      

      //vector from ith to jth wannier center
      D3vector r0 = moli -molj;
      //D3vector r0 = D3vector (5,5,5);;
      c_.images(r0);
      double w[9];

      //include interaction of neighboring cells
      int ncell = 40;
      double eps = 0.001;
      int ncount = 0;
      for ( int ix = - ncell; ix <= ncell; ix ++)
      for ( int iy = - ncell; iy <= ncell; iy ++)
      for ( int iz = - ncell; iz <= ncell; iz ++)
      { 
  
        D3vector vec_cell ( c_.x() * ix, iy * c_.y(), iz * c_.z() );

        D3vector r = r0 + vec_cell;


        double rr = length (r);     
        if ( rr < eps ) continue;

        double r3 = rr * rr * rr;
        double r5 = r3 * rr * rr;

        if ( rr - eps >  c_.x() * ( ncell + 0.5 ) ) continue;
        //construct interaction matrix W
        for ( int ii = 0; ii < 3; ii ++)
          for ( int jj = 0; jj < 3; jj ++)
        {
          int ind = ii * 3 + jj;

          w [ ind ] = 3.0 * r[ii] * r[jj] / r5;

          if ( ii == jj ) w[ind] -= 1.0 / r3;
        }
        ncount ++;
#if 0
        cout << r  << "  " << rr << endl;
        for ( int ii = 0; ii < 3; ii ++ )
          for ( int jj = 0; jj < 3; jj ++ )
        //if ( jj >= ii )
          cout << setw(12) << w[jj+ii*3] 
               << setw(12) << polarj[jj+ii*3] << endl;
#endif

        //elocal += W * p = efield_ * W * alpha_eff
        dgemm( &cn, &cn, &three, &three, &three, &one, w, 
                &three, polarj, &three, &one, elocal_loc, &three);
      }
    
      #pragma omp critical
      for ( int ii = 0; ii < 9; ii ++ ) elocal[ii] += elocal_loc[ii];

      }//end pragma omp parallel


    }// for j

    //elocal += 4pi / 3  * P
    //daxpy ( & nine, &depolar, polari, &ione, elocal, &ione );

    //local field
    daxpy ( & nine, &one, eext_, &ione, elocal, &ione );   

    //copy w to w_inverse and inverse
    elocal_inv_.init(elocal,3);
    elocal_inv_.inverse();

    //alpha = dp * W^-1    
    dgemm( &cn, &cn, &three, &three, &three, &one, polari, &three, 
          elocal_inv_.valptr(), &three, &zero, molset_[i]->alpha_tensor(), &three);
  
    //cout << "Calc done\n";

#if 0
         
    //test
    double test[9];
    double zero = 0.0;
    dgemm( &cn, &cn, &three, &three, &three, &one, molset_[i].alpha_tensor(), 
      &three, elocal, &three, &zero, test, &three);
  
    if ( i == 0 )
    for ( int ii = 0; ii < 9; ii ++ ) 
      cout << setw(15) << eext_[ii] 
           << setw(15) << elocal[ii]
           << setw(15) << elocal[ii] - eext_[ii] 
           << setw(15) << polari[ii]
           << setw(15) << elocal_inv_[ii]
           << setw(15) << molset_[i].alpha_tensor()[ii] 
           << setw(15) << polari[ii] 
           << setw(15) << test[ii] << endl;
#endif

  }//for i


}


void Polariz::Calculate_mlwf(void)
{

  for ( int i = 0; i < 1; i ++ )//mod
  {

    D3vector wfi = mlwfset_[i].x();

    //effective polarizability tensor
    double * polari = mlwfset_[i].polar_tensor();
    
    int nine = 9, ione = 1, three = 3; 
    double one = 1.0;
    char cn = 'n';

    //elocal restore local field of WF i
    double  elocal[9];
    for ( int ii = 0; ii < 9; ii ++ ) elocal[ii] = 0.0;
  
    for ( int j = 0; j < nmlwf_; j ++ ) //mod 
    {
      if ( i == j ) continue;

      const D3vector& wfj = mlwfset_[j].x();

      //effective polarizability tensor
      double * polarj = mlwfset_[j].polar_tensor();

      //vector from ith to jth wannier center
      D3vector r0 = wfi -wfj;
      c_.images(r0);
      double w[9];

      //include interaction of neighboring cells
      int ncell = 5;
      for ( int ix = - ncell; ix <= ncell; ix ++)
      for ( int iy = - ncell; iy <= ncell; iy ++)
      for ( int iz = - ncell; iz <= ncell; iz ++)
      { 
        D3vector vec_cell ( c_.x() * ix, iy * c_.y(), iz * c_.z() );

        D3vector r = r0 + vec_cell;
        double rr = length (r);     

        if ( rr >  c_.x() * ( .5 + ncell ) ) continue;
        //construct interaction matrix W
        for ( int ii = 0; ii < 3; ii ++)
          for ( int jj = 0; jj < 3; jj ++)
        {
          int ind = ii * 3 + jj;

          w [ ind ] = -3.0 * r[ii] * r[jj] / pow ( rr, 5 );

          if ( ii == jj ) w[ind] += 1.0 / pow ( rr, 3 );
        }
#if 1
        cout << r0  << "  " << rr << "  ";
        for ( int ii = 0; ii < 1; ii ++ )
          for ( int jj = 0; jj < 3; jj ++ )
        if ( jj >= ii )
          cout << setprecision(8) << setw(10) << w[jj] << "  ";
        cout << endl;
#endif

        //elocal += W * p
        dgemm( &cn, &cn, &three, &three, &three, &one, w, 
                &three, polarj, &three, &one, elocal, &three);
      }

    }// for j

    //local field
    daxpy ( & nine, &one, eext_, &ione, elocal, &ione );   

    //copy w to w_inverse and inverse
    Context ctxt;
    DoubleMatrix e_inverse( ctxt, 3, 3);
    e_inverse.init(elocal,3);
    e_inverse.inverse();

    //alpha = dp * W^-1    
    dgemm( &cn, &cn, &three, &three, &three, &one, polari, &three, 
          e_inverse.valptr(), &three, &one, mlwfset_[i].alpha_tensor(), &three);

    for ( int ii = 0; ii < 9; ii ++ ) cout << mlwfset_[i].alpha_tensor()[ii] << endl;
    cout << endl;

  }//for i

}


Polariz::~Polariz()
{
  delete ewald_mol_;
  delete ewald_mlwf_;
}
