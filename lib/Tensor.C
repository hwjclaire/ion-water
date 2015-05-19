////////////////////////////////////////////////////////////////////////////////
//
// Tensor.h
//
// double 3-vectors
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "blas.h"
#include "D3vector.h"
#include "Tensor.h"

using namespace std;

Tensor::Tensor(double* rhs)
{
  for ( int i = 0; i < 9; i ++ )
    a_[i] = rhs[i];
}


Tensor::Tensor (const double& xx, const double& yy, const double& zz,
                const double& xy, const double& yz, const double& zx,
                const char & uplo)
{ 
  a_[0] = xx; 
  a_[4] = yy;
  a_[8] = zz;

  if ( uplo == 'l' )
  {
    a_[1] = xy;
    a_[2] = zx;
    a_[5] = yz;
  }
  else if ( uplo == 'u' )
  {
    a_[3] = xy;
    a_[6] = zx;
    a_[7] = yz;
  }
  else if ( uplo == 's' )
  {
    a_[1] = xy;
    a_[2] = zx;
    a_[5] = yz;
    a_[3] = xy;
    a_[6] = zx;
    a_[7] = yz;
  }
  else 
    assert(false);

}
  
Tensor::Tensor(const D3vector& diag, const D3vector& offdiag)
{
  a_[0] = diag[0];
  a_[4] = diag[1];
  a_[8] = diag[2];
  a_[1] = offdiag[0];
  a_[5] = offdiag[1];
  a_[2] = offdiag[2];
  a_[3] = offdiag[0];
  a_[7] = offdiag[1];
  a_[6] = offdiag[2];
}

Tensor::Tensor(const double * diag, const double * offdiag)
{
  a_[0] = diag[0];
  a_[4] = diag[1];
  a_[8] = diag[2];
  a_[1] = offdiag[0];
  a_[5] = offdiag[1];
  a_[2] = offdiag[2];
  a_[3] = offdiag[0];
  a_[7] = offdiag[1];
  a_[6] = offdiag[2];
}



void Tensor::setdiag(const int& i, const double& b)
{
  assert(i>=0&&i<3);
  a_[i*4] = b;
}

void Tensor::setdiag(const D3vector& b)
{
  for ( int i = 0; i < 3; i ++ )
    a_[i*4] = b[i];
}

void Tensor::setoffdiag(const int& i, const double& b)
{
  assert(i>=0&&i<3);
  if ( i == 0 )
  {
    a_[1] = b;
    a_[3] = b;
  }
  else if ( i == 2 )
  {
    a_[2] = b;
    a_[6] = b;
  }
  else
  {
    a_[5] = b;
    a_[7] = b;
  }
}

void Tensor::setoffdiag(const D3vector& b)
{
  a_[1] = b[0];
  a_[3] = b[0];
  a_[5] = b[1];
  a_[7] = b[1];
  a_[2] = b[2];
  a_[6] = b[2];
}

Tensor operator * ( Tensor& a, Tensor& b)
{
  Tensor c;
  int ithree = 3, ione = 1;
  double one = 1.0, zero = 0.0;
  char t = 'n';
  dgemv ( &t, &ithree, &ithree, &one, &a[0], &ithree,
          &b[0], &ione, &zero, &c[0], &ione );
  return c;
}

D3vector operator * ( Tensor& a, D3vector& b)
{
  D3vector c;
  int ithree = 3;
  double one = 1.0, zero = 0.0;
  char t = 'n';
  dgemm ( &t, &t, &ithree, &ithree, &ithree, &one, &a[0], &ithree,
          &b[0], &ithree, &zero, &c[0], &ithree );
  return c;
}



double Tensor::trace(void)
{
  return a_[0]+a_[4]+a_[8];
}

void Tensor::traceless(void)
{
  double b = trace() / 3;
  a_[0] -= b;
  a_[4] -= b;
  a_[8] -= b;
}

void Tensor::clear(void)
{
  for ( int i = 0; i < 9; i ++ )
    a_[i] = 0.0;
}

void Tensor::identity(void)
{
  clear();
  a_[0] = 1.0;
  a_[4] = 1.0;
  a_[8] = 1.0;
}

std::ostream& operator << ( std::ostream& os, const Tensor& rhs )
{
  const double * const v  = rhs.a();
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os.precision(6);
  os << setw(12) << v[0] << setw(12) << v[3] << setw(12) << v[6] << "\n"
     << setw(12) << v[1] << setw(12) << v[4] << setw(12) << v[7] << "\n"
     << setw(12) << v[2] << setw(12) << v[5] << setw(12) << v[8] << "\n";
  return os;
}

void Tensor::print_inline( std::ostream& os )
{
  const double * const v  = &a_[0];
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os.precision(10);
  os << setw(16) << v[0] << setw(16) << v[4] << setw(16) << v[8] 
     << setw(16) << v[1] << setw(16) << v[5] << setw(16) << v[2];

}


void Tensor::syev(char uplo, D3vector& eigval, Tensor& eigvec)
{
  double w[3];

  int info;
  char jobz = 'V';
  int lwork=-1;
  double tmplwork;
  int n = 3;

  dsyev(&jobz, &uplo, &n, eigvec.a(), &n, &w[0], &tmplwork, &lwork, &info);

  lwork = (int) tmplwork + 1;
  double* work = new double[lwork];

  eigvec = *this;
  dsyev(&jobz, &uplo, &n, eigvec.a(), &n, &w[0], work, &lwork, &info);
  delete[] work;

  eigval = D3vector(&w[0]);

}
