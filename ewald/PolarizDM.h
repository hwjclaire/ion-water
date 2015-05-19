#ifndef POLARIZ_H
#define POLARIZ_H

#include <vector>
#include <iostream>
#include <cassert>
#include "Context.h"
#include "Matrix.h"
#include "Mlwf.h"
#include "Cell.h"
#include "D3vector.h"
#include "Ewald.h"

using namespace std;


class Polariz
{
  private:

  Cell& c_; //Unit cell

  Context ctxt_;
  DoubleMatrix elocal_inv_;

  std::vector < Mlwf > & mlwfset_;


  std::vector < D3vector * >  mlwf_pos_;
  std::vector < double * >  mlwf_spread_;
  std::vector < double * >  mlwf_polar_;
  std::vector < std::vector < double > >  mlwf_elocal_;

  const int nmol_;
  const int nmlwf_;

  double eext_[9];

  //for ewald:
  Ewald *ewald_mol_, *ewald_mlwf_ ;
  double kappa_;
  int nk_;

  public:

  Polariz(Cell & c, std::vector < Mol* > & molset, std::vector < Mlwf > & mlwfset_);

  void Calculate_mol_Ewald(void);

  void Calculate_mol(void);

  void Calculate_mlwf_Ewald(void);

  void Polar_mlwf_rotate(void);

  void Calculate_mlwf(void);

  ~Polariz();

};

#endif 
