#ifndef STAT_H
#define STAT_H

#include <iostream>
#include <cassert>
#include <vector>

using namespace std;

class Stat
{

  private:

  double  sum_, sum2_, max_, min_, bin_;
  int n_;

  double vmax_, vmin_;//max & min values

  int nbin_ ; 

  std::vector<double> distr_;

  public:

  Stat(double min, double max);

  Stat(double min, double max, double bin_);

  Stat(const Stat& s);

  void init(double min, double max, double bin);

  void init(double min, double max);

  double avg() const { assert (n_ > 0); return sum_ / n_; }

  double rmin() const {return min_; }

  double rmax() const {return max_; }

  double sum() const {return sum_; }

  double sum2() const { return sum2_; }
 
  double max() const { return vmax_;}
  
  double min() const { return vmin_;}

  double bin() const { return bin_;}

  int n() const { return n_; }

  double stddev() const;

  void reset();

  void add(double x);

  void print(ostream& os);
};
#endif
