#ifndef NG_H
#define NG_H

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
class NG: public Mol
{

  private:
  
  int natom_;
  
  Atom *o_ ;

  void add_sulfur( Atom& sulfur) {};

  void add_hydrogen( Atom& hydrogen) {};

  void add_halide( Atom& hydrogen) {};

  public:
 
  Atom* o() { return o_; } 
  int& number() { return number_; }
  bool wf_full() { return nwf_ == 4; }
  bool atom_full(){ return natom_ == 1 ; }
  

  NG(int n , Cell& c)
  {
    number_ = n;
    c_ = &c;
    wf_.resize(4,0);
  }

  void reset_wf(void);

  void reset_atom();

  void add_wf( Mlwf & a);

  void add_oxygen( Atom & oxygen);

  void check_atom ();
  
  void check_wf();

  double distance( Atom * );
  
  double distance( Mlwf * );

  void Compute_cm(); 
   
  void Compute_cc(); 
   
  void Compute_dipole();
};
#endif
