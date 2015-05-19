////////////////////////////////////////////////////////////////////////////////
//
// Tensor.h
//
// double 3-vectors
//
////////////////////////////////////////////////////////////////////////////////

#ifndef TENSOR_H
#define TENSOR_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "blas.h"
#include "D3vector.h"

using namespace std;

class Tensor
{
  private:

  double a_[9];

  public:
  
  double* a (void) { return &a_[0];}
  const double* a (void) const { return &a_[0];} 

  explicit Tensor(double* rhs);

  explicit Tensor(void) {clear();}

  explicit Tensor(const double& xx, const double& yy, const double& zz)
  { a_[0]=xx; a_[4]=yy; a_[8]=zz; }

  explicit Tensor(const double& xx, const double& yy, const double& zz,
                  const double& xy, const double& yz, const double& zx,
                  const char & uplo);
  
  explicit Tensor(const D3vector& diag, const D3vector& offdiag);

  explicit Tensor(const double * diag, const double * offdiag);

  explicit Tensor(const double * a)
  {
    for ( int i = 0; i < 9; i ++ ) a_[i] = a[i];
  }

  double& operator[](const int &i)
  {
    assert(i>=0 && i < 9);
    return a_[i];
  }

  const double& operator[](const int &i) const 
  {
    assert(i>=0 && i < 9);
    return a_[i];
  }

  void setdiag(const int& i, const double& b);

  void setdiag(const D3vector& b);

  void setoffdiag(const int& i, const double& b);

  void setoffdiag(const D3vector& b);

  bool operator==(const Tensor &rhs) const
  {
    bool eq = true;
    for ( int i = 0; i < 9; i ++ )
    {
      if ( rhs[i] == a_[i] )
      {
        eq = false;
        break;
      }
    }
    return eq;
  }

  bool operator!=(const Tensor &rhs) const
  {
    bool neq = false;
    for ( int i = 0; i < 9; i ++ )
    { 
      if ( rhs[i] != a_[i] )
      {
        neq = true;
        break;
      }
    }
    return neq;
  }

  Tensor& operator += ( const Tensor& rhs )
  {
    for ( int i = 0; i < 9; i ++ )
      a_[i] += rhs[i];
    return *this;
  }

  Tensor& operator -= ( const Tensor& rhs )
  {
    for ( int i = 0; i < 9; i ++ )
      a_[i] -= rhs[i];
    return *this;
  }

  Tensor& operator *= ( const double& rhs )
  {
    for ( int i = 0; i < 9; i ++ )
      a_[i] *= rhs;
    return *this;
  }

  Tensor& operator /= ( const double& rhs )
  {
    for ( int i = 0; i < 9; i ++ )
      a_[i] /= rhs;
    return *this;
  }

  friend const Tensor operator + (const Tensor& lhs, const Tensor& rhs )
  {
    return Tensor(lhs) += rhs;
  }

  friend const Tensor operator - ( const Tensor& a, const Tensor& b )
  {
    return Tensor(a) -= b;
  }

  friend Tensor operator * ( const double& a, const Tensor& b )
  {
    return Tensor(b) *= a;
  }

  friend Tensor operator * ( const Tensor& a, const double& b )
  {
    return Tensor(a) *= b;
  }

  friend Tensor operator / ( const Tensor& a, const double& b )
  {
    return Tensor(a) /= b;
  }

  friend Tensor operator * ( Tensor& a, Tensor& b);

  friend D3vector operator * ( Tensor& a, D3vector& b);

  friend Tensor operator - ( const Tensor& a ) // unary minus
  {
    return Tensor()-a;
  }

  double norm2( const Tensor& a )
  {
    int inine = 9, ione = 1;
    return dnrm2( &inine, &a[0], &ione );
  }

  double norm( const Tensor& a )
  {
    return sqrt ( norm2 (a) );
  }

  double trace(void);
  
  void traceless(void);

  void clear(void);

  void identity(void);

  friend std::ostream& operator << ( std::ostream& os, const Tensor& rhs );

  void print_inline( std::ostream& os );

  void syev(char uplo, D3vector& eigval, Tensor& eigvec);

};
#endif
