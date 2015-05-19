#ifndef WATER_H
#define WATER_H

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Atom.h"
#include "Mlwf.h"
#include "D3vector.h"
#include "Mol.h"

class Water: public Mol
{

  private:
  
  Atom *o_, *h_[3] ;

  int nh_;

  void add_sulfur( Atom& sulfur) {};

  void add_halide( Atom& sulfur) {};

  void add_alkali( Atom & alkali){};

  //primary axes:
  //1. dipole axis
  //2. axis normal to molecule plain
  //3. axis parallel to molecule plain
  D3vector axis_[3];
  double axis_tensor_[9]; 

  D3vector oh_[2]; // vector that point to hydrogen from oxygen 
  D3vector owf_[4]; // vector point to WFc from oxygen

  double alpha_tensor_axis_[9];  
  double polar_tensor_axis_[9];  

  public:
  
  Atom* o() const { return o_;}
  Atom* h(int i) const { return h_[i];}
  int& number() { return number_; }
  bool isoh(){ return nh_ == 1; }
  bool isoh2(){ return nh_ == 2; }
  bool isoh3(){ return nh_ == 3; }
  bool wf_full() { return nwf_ == 4; }
  bool atom_full(){ return h_[2] != 0; }
  D3vector& axis_vec(int n) { return axis_[n];}
  D3vector& owf(int n) { assert(n>=0&&n<nwf_); return owf_[n];}
  D3vector& oh(int n) { assert(n>=0&&n<nh_); return oh_[n];}
  double* axis_ptr(int n) { return &axis_tensor_[3*n];}
  double* axis_tensor() { return &axis_tensor_[0];}
  double* alpha_tensor_axis() {return alpha_tensor_axis_;}
  double& alpha_tensor_axis(int n) {return alpha_tensor_axis_[n];}
  double* polar_tensor_axis() {return polar_tensor_axis_;}
  double& polar_tensor_axis(int n) {return polar_tensor_axis_[n];}
  double * alpha_tensor(void) {return &alpha_tensor_[0];}
  double * polar_tensor(void) {return &polar_tensor_[0];}



  Water(int n , Cell& c) : o_(0),nh_(0)
  {
    number_ = n;
    c_ = &c;
    for ( int i = 0; i < 3; i ++) h_[i] = 0;
    wf_.resize(4,0);
  }

  void reset_wf(void);

  void reset_atom();

  void add_wf( Mlwf & a, bool noorder);

  void add_wf( Mlwf & a) { add_wf(a,true); };

  void add_oxygen( Atom & oxygen);

  void add_hydrogen( Atom & hydrogen );

  void check_atom ();
  
  void check_wf();

  void getaxis();

  void getaxis_dipole();

  void rotate_polar_tensor();
  
  void rotate_alpha_tensor();

  void watermove();

  void mlwfmove();

  double distance( Atom * );
  
  double distance( Mlwf * );

  void Compute_cm(); 
   
  void Compute_cc(); 
   
  void Compute_dipole();

  void Compute_quad();

  void Compute_polariz();
};
#endif
