#ifndef POLARIZ_H
#define POLARIZ_H

#include <vector>
#include <iostream>
#include <cassert>
#include "Context.h"
#include "Mlwf.h"
#include "Cell.h"
#include "D3vector.h"
#include "Ewald.h"
#include "Tensor.h"

using namespace std;


class Polariz
{
  private:

  Cell& c_; //Unit cell

  std::vector < Tensor > dm_;

  std::vector < D3vector >  mlwf_pos_;
  std::vector < double >  mlwf_polar_;

  const int nmlwf_;

  //for ewald:
  Ewald *ewald_;
  double kappa_;
  int nk_;

  public:

  Polariz(Cell & c, std::vector < Mlwf > & mlwfset);

  const std::vector < Tensor > & dm(void) const {return dm_;};

  void Calculate_dm(void);

  ~Polariz();

};

#endif 
