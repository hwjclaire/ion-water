#ifndef MOL_H
#define MOL_H
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include "Cell.h"
#include "Atom.h"
#include "Mlwf.h"
#include "D3vector.h"
#include "Tensor.h"

class Mlwf;
class Atom;

class Mol
{

  private:

  protected:

  std::vector < Mlwf* > wf_; //wannier centers

  int nwf_; //number of wfs

  int number_; //numbering in a set of molecules

  Atom * center_; //pointer to the central atom, sulfur for sulfate and oxygen for water

  Cell * c_; // pointer to a unitcell
  
  D3vector dipole_ , cm_, cc_ ;//dipole and cneter of mass and center of charge

  double tot_charge_;//total charge

  double alpha_tensor_[9];
  double polar_tensor_[9];

  Tensor polariz_;
  Tensor quad_, quad_t_, quad_ion_, quad_mlwfc_, quad_mlwfs_;
  D3vector quad_eigval_;
  Tensor quad_eigvec_;  

  public:
  
  double * alpha_tensor(void) {return &alpha_tensor_[0];}
  double * polar_tensor(void) {return &polar_tensor_[0];}

  int& nwf() { return nwf_; }

  int& number() { return number_; }

  Mlwf* wf(int n)
  {
    assert ( n < nwf_ );
    return wf_[n];
  }
  
  Mol() : nwf_(0) {};

  double tot_charge() { return tot_charge_;}

  D3vector& dipole() { return dipole_; }
  D3vector& cm() { return cm_; }
  D3vector& cc() { return cc_; }

  Tensor& polariz() { return polariz_; }
  Tensor& quad() { return quad_; }
  Tensor& quad_t() { return quad_t_; }
  D3vector& quad_eigval() { return quad_eigval_; }
  Tensor& quad_eigvec() { return quad_eigvec_; }
  Tensor& quad_ion() { return quad_ion_; }
  Tensor& quad_mlwfc() { return quad_mlwfc_; }
  Tensor& quad_mlwfs() { return quad_mlwfs_; }

  virtual void reset_atom(void) = 0;

  virtual void reset_wf(void) = 0;

  virtual void add_wf( Mlwf & a, bool noorder) = 0;

  virtual void add_wf( Mlwf & a) = 0;

  virtual void add_sulfur( Atom & oxygen) = 0;

  virtual void add_oxygen( Atom & oxygen) = 0;

  virtual void add_hydrogen( Atom & hydrogen ) = 0;

  virtual void add_halide( Atom & halide ) = 0;

  virtual void add_alkali( Atom & ) = 0;

  virtual void check_atom () = 0;
  
  virtual void check_wf() = 0;

  virtual bool atom_full() = 0;

  virtual bool wf_full() = 0;

  virtual double distance( Mlwf * ) = 0;
  
  virtual double distance( Atom * ) = 0;

  virtual void Compute_dipole() = 0;

  virtual void Compute_quad() = 0;

  virtual void Compute_cm() = 0;

  virtual void Compute_cc() = 0;

  virtual void Compute_polariz() = 0;

};
#endif
