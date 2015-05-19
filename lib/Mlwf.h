#ifndef MLWF_H
#define MLWF_H
 
#include "D3vector.h"
#include "iostream"
#include "Mol.h"
#include "Tensor.h"

class Mol;
class Mlwf
{
  private:

  int number_;//number in a moleucle.

  Mol * mol_ ;

  D3vector x_;
  
  double spread_;

  double polar_[6];

  double polar_tensor_[9];

  double alpha_tensor_[9];

  double axis_tensor_[9];

  double alpha_tensor_axis_[9];

  double polar_tensor_axis_[9];

  Tensor quad_;
  Tensor polariz_;

  public:

  
  void setnumber(int n ) { number_ = n;};
  int number() { return number_;};
  Mlwf(void) { mol_ = 0; }
  Mol* mol(void) { return mol_; }
  void setmol( Mol* rhs ) { mol_ = rhs; }
  double * polar() { return &polar_[0] ;}; 
  double * alpha_tensor() { return &alpha_tensor_[0];}; 
  double * polar_tensor() { return &polar_tensor_[0];}; 
  double * axis_tensor() { return &axis_tensor_[0];}; 
  double * alpha_tensor_axis() { return &alpha_tensor_axis_[0];}; 
  double * polar_tensor_axis() { return &polar_tensor_axis_[0];}; 


  D3vector& x() {return x_;}
  double& spread() {return spread_;}
  double spread2() {return spread_*spread_;}

  Tensor& quad() { return quad_; }  
  Tensor& polariz() { return polariz_; }  

  void readx(istream& is) {is >> x_ ;}

  void readpos(istream& is) {is >> x_;}
  
  void reads(istream& is) {is >> spread_;}
  
  void readp(istream& is)
  {
    for ( int i = 0; i < 6; i++) is >> polar_[i];
  
    polar_tensor_[0]=polar_[0];    
    polar_tensor_[4]=polar_[3];    
    polar_tensor_[8]=polar_[5];    
    polar_tensor_[1]=polar_[1];    
    polar_tensor_[5]=polar_[2];    
    polar_tensor_[2]=polar_[4];     
    polar_tensor_[3]= polar_tensor_[1];
    polar_tensor_[6]= polar_tensor_[2];
    polar_tensor_[7]= polar_tensor_[5];

    polariz_ = Tensor(polar_tensor_);
  }

  void readquad(istream& is)
  {
    double xx,yy,zz,xy,yz,zx;

    is >> xx >> yy >> zz >> xy >> yz >> zx;

    quad_ = Tensor(xx,yy,zz,xy,yz,zx,'s');

  }

};
#endif
