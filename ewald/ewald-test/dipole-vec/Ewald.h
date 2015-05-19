
#ifndef EWALD_H
#define EWALD_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
#include "D3vector.h"
#include "Cell.h"


using namespace std;

class Ewald
{
  private:

  Cell& cell_;

  //constants: std dev of cancelling Gaussian charge
  double kappa_, kappa2_, kappa3_, fourkappa2_, const_self_;
  const double pi_;
  const double sqrtpi_;

  int nkx_, nky_, nkz_, nk_; //number of k vectors in each dimension

  std::vector < D3vector > k_; // k vectors
  
  std::vector < double > k2_; // norm2 of k vectors

  std::vector < std::complex < double > > rho_k_; // charge density in reciprocal space

  int nion_;

  std::vector < D3vector > ion_; // position of ions

  std::vector < std::vector < double > > ionpolar_; // polarization of ions

  std::vector < std::vector < double > >  ionfield_short_; // short range potential of ions

  std::vector < std::vector < double > > ionfield_long_; // long range potential of ions

  std::vector < std::vector < double > > ionfield_self_; // self potential of ions

  public:
  
  Ewald(Cell& cell, int nk, double kappa);

  //Ewald(int nion, std::vector<Mol*> mol);

  void potential(void){};

  void shortrange(void); // short range
  
  void longrange(void); // long range
  void longrange1(void); // long range

  void self(void); // self potential

  void print ( std::ostream& os);

};
#endif
