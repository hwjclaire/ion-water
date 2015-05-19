
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

  double kappa_, fourkappa2_;//std dev of cancelling Gaussian charge

  int nkx_, nky_, nkz_, nk_; //number of k vectors in each dimension

  std::vector < D3vector > k_; // k vectors
  
  std::vector < double > k2_; // norm2 of k vectors

  std::vector < std::complex < double > > sk_; // long range potential in reciprocal space

  int nion_;

  std::vector < D3vector > ion_; // position of ions

  std::vector < double > ioncharge_; // charge of ions

  std::vector < double > ionpotential_short_; // short range potential of ions

  std::vector < double > ionpotential_long_; // long range potential of ions

  std::vector < double > ionpotential_self_; // self potential of ions

  public:
  
  Ewald(Cell& cell);

  //Ewald(int nion, std::vector<Mol*> mol);

  void potential(void){};

  void shortrange(void); // short range
  
  void longrange(void); // long range
  void longrange1(void); // long range

  void self(void); // self potential

  void print ( std::ostream& os);

};
#endif
