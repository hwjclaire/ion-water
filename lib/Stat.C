#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "Stat.h"

using namespace std;




Stat::Stat(double min, double max)
{
  init(min,max);
}

Stat::Stat(double min, double max, double bin)
{
  init(min,max,bin);
}

Stat::Stat(const Stat& s)
{
  init(s.rmin(),s.rmax(),s.bin());
}

void Stat::init(double min, double max)
{
  init(min, max, 0.01);
}

void Stat::init(double min, double max, double bin)
{
  sum_=0.0;
  sum2_=0.0; 
  n_=0;
  max_=max;
  min_=min;
  bin_=bin;

  nbin_ = ( max_ - min_ ) / bin_;

  assert ( max_ > min_ );
  assert ( nbin_ > 0 );

  distr_.resize(nbin_,0.0);

  //redefine max_
  max_ = min_ + bin_ * nbin_;
}

double Stat::stddev() const
{ 
  assert ( n_ > 0 );
  return sqrt ( sum2_ / n_ - sum_ * sum_ / n_ / n_ );
}

void Stat::reset()
{
  n_ = 0; sum_ = 0; sum2_ = 0;  
  for ( int i = 0; i < nbin_; i++) distr_[i] =0;
}  

void Stat::add(double x)
{
  if ( x> max_ || x < min_ ) cout << "Stat out of range: " << x << endl;

  //assert ( x > min_ && x < max_ );

  if ( n_ == 0 )
  {
    vmax_ = x; 
    vmin_ = x;
  }
  else
  {
    if ( x > vmax_ ) vmax_ = x;
    if ( x < vmin_ ) vmin_ = x;
  }

  sum_ += x; sum2_ += x*x;

  for ( int i = 0; i < nbin_; i++)
  {
    if ( x < i * bin_ + min_ ) 
    {
      distr_[i-1] += 1.0;
      break;
    }
  }

  n_++;    

}

void Stat::print(ostream& os)
{
  os.precision(6);
  for ( int i = 0; i < nbin_; i++)
  {
    os << setw(12) << min_ + bin_ * ( i + 0.5 ) << setw(12) << distr_[i] / n_ / bin_ << endl;
  }
}

