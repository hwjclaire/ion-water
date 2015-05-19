#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include "D3vector.h"
#include "Cell.h"
#include "Mlwf.h"
#include "Mol.h"

class Mlwf;
class Mol;
class Atom
{
  private:

  double mass_, charge_;

  D3vector x_,v_;

  Mol * mol_;
  
  string name_;

  public:

  Atom(void) { mol_ = 0; }
  double charge() { return charge_;}
  void setcharge(double rhs) { charge_ = rhs;}
  double mass(void) { return mass_;}
  void setmass(double rhs) { mass_ = rhs; }
  Mol* mol(void) { return mol_; }
  void setmol( Mol * const rhs ) { mol_ = rhs; }
  string& name() { return name_;}
  void setname(const string& rhs) { name_ = rhs;}
  //void setname(char* rhs) { name_ = rhs;}

  void reset(){mol_ = 0;}

  const D3vector& x(void) const { return x_; };
  const D3vector& v(void) const { return v_; };

  friend std::istream& operator >> ( std::istream& is, Atom& a )
  {
    is >> a.x_ >> a.v_ ;
    return is;
  }

  void readx ( std::istream& is ) { is >> x_; }
  
  void readv ( std::istream& is ) { is >> v_; }

  double distance( Atom *a, const Cell * const c)
  {
    D3vector vec = a -> x() -x_;
    c -> images(vec);
    return length(vec);
  }
  
  double distance( Mlwf *a, const Cell * const c)
  {
    D3vector vec = a -> x() -x_;
    c -> images(vec);
    return length(vec);
  }


};
#endif
