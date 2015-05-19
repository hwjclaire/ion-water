#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Atom.h"
#include "Water.h"
#include "Mlwf.h"
#include "D3vector.h"
#include "Mol.h"

class Sulfate : public Mol
{

  private:
  
  Atom *s_, *o_[4], *h1_, *h2_;

  bool minus2_; bool minus1_;

  void add_halide( Atom& hydrogen) {};

  void add_alkali( Atom & alkali){};

  public:
  
  Atom* s() { return s_;}
  Atom* o(int n) { return o_[n];}
  Atom* h1() { return h1_;}
  Atom* h2() { return h2_;}
  int& number() { return number_; }
  bool minus1(){ return h2_ == 0 && h1_ != 0; }
  bool minus2(){ return h2_ == 0 && h1_ == 0; }
  bool wf_full() { return nwf_ == 16; }
  bool atom_full() { return h2_ != 0; }
  const D3vector& cm() { return cm_; }

  Sulfate(int n,Cell& c) : s_(0), h1_(0), h2_(0),
     minus1_(false), minus2_(false)
  {
    number_ = n;
    c_ = &c;
    for ( int i = 0; i < 4; i++) o_[i] = 0;
    wf_.resize(16,0);
  }

  void reset_atom(void);

  void reset_wf(void);

  void add_wf( Mlwf & a, bool noorder);

  void add_wf( Mlwf & a);

  void add_sulfur( Atom & oxygen);

  void add_oxygen( Atom & oxygen);

  void add_hydrogen( Atom & hydrogen );

  void check_atom ();
  
  void check_wf();

  double distance( Atom * ) ;

  double distance( Mlwf * ) ;

  double distance( Water * a ) ;
  
  void Compute_cm();

  void Compute_cc(); 

  void Compute_dipole();

  void Compute_quad();

  void Compute_polariz();
};
