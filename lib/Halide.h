#ifndef HALIDE_H
#define HALIDE_H

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Atom.h"
#include "Mlwf.h"
#include "D3vector.h"
#include "Mol.h"

//noble gas
class Halide: public Mol
{

  private:
  
  int natom_;
  
  Atom *xx_ ;

  void add_sulfur( Atom& sulfur) {};

  void add_hydrogen( Atom& hydrogen) {};

  void add_oxygen( Atom & oxygen){};

  void add_alkali( Atom & alkali){};

  void Compute_quad(){};

  public:
 
  Atom* xx() { return xx_; } 
  int& number() { return number_; }
  bool wf_full() { return nwf_ == 4; }
  bool atom_full(){ return natom_ == 1 ; }

  Halide(int n , Cell& c)
  {
    number_ = n;
    c_ = &c;
    wf_.resize(4,0);
  }

  void reset_wf(void);

  void reset_atom();

  void add_wf( Mlwf & a, bool);

  void add_wf( Mlwf & a) { add_wf(a,true); };

  void add_halide( Atom& hydrogen);

  void check_atom ();
  
  void check_wf();

  double distance( Atom * );
  
  double distance( Mlwf * );

  void Compute_cm(); 
   
  void Compute_cc(); 
   
  void Compute_dipole();

  void Compute_polariz();
};
#endif
