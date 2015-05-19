#ifndef CELL_H
#define CELL_H

#include "D3vector.h"
#include <cassert>
#include <vector>


using namespace std;

class Cell
{
  private:
  
  std::vector < double > a_, b_; 
  double volume_;

  public:
 
  double x(){return a_[0];}
  double y(){return a_[4];}
  double z(){return a_[8];}
  
  double v(){return volume_;} 

  const std::vector < double >& a(void) const {return a_;}

  Cell(void) {a_.resize(9);b_.resize(9);}

  Cell ( const std::vector < double > & a ) 
  { 
    b_.resize(9);
    assert ( a.size() == 9 );
    a_ = a; 
    invert();
  }

  friend std::istream& operator >> ( std::istream& is, Cell& c )
  {
    for ( int i = 0; i < 9; i ++)
      is >> c.a_[i];
    return is;
  }

  friend std::ostream& operator << ( std::ostream& os, Cell& c )
  {
    os << "cell:" << endl;
    for ( int i = 0; i < 9; i ++)
      os << c.a_[i] << endl;

    os << "rcell:" << endl;
    for ( int i = 0; i < 9; i ++)
      os << c.b_[i] << endl;

    return os;
  }


  void invert(void)
  {
    assert ( a_.size() == 9);
    assert ( b_.size() == 9);

    //adjoint matrix
    b_[0]=a_[4]*a_[8]-a_[5]*a_[7];
    b_[3]=a_[7]*a_[2]-a_[8]*a_[1];
    b_[6]=a_[1]*a_[5]-a_[2]*a_[4];
    b_[1]=a_[5]*a_[6]-a_[8]*a_[3];
    b_[4]=a_[8]*a_[0]-a_[2]*a_[6];
    b_[7]=a_[2]*a_[3]-a_[5]*a_[0];
    b_[2]=a_[3]*a_[7]-a_[6]*a_[4];
    b_[5]=a_[6]*a_[1]-a_[0]*a_[7];
    b_[8]=a_[0]*a_[4]-a_[3]*a_[1];
 
    //determinant

    volume_ = a_[0]*b_[0]+a_[3]*b_[3]+a_[6]*b_[6];

    //
    for(int i = 0; i < 9; i ++)
    { 
      b_[i] = b_[i] / volume_;
    }

  } //invert



  void images( D3vector& vec) const
  {

    double t[3];

    for(int i = 0; i < 3; i++)
    {
      t[i]=0.;
      for( int j = 0; j < 3; j++) 
          t[i]+= vec[j] * b_ [ j + 3 * i ];
    
      t[i]-=(int)t[i];
      /////////////////fold to from -0.5a to 0.5a
      if(t[i]>0.5)t[i]=t[i]-1;
      if(t[i]<-0.5)t[i]=t[i]+1;
    }
 
    for(int i = 0; i < 3; i++)
    {
      vec[i]=0;
      for(int j = 0; j < 3; j++) 
          vec[i] += t[j] * a_ [ j + 3 * i ];
    }

  }//images

};//cell

#endif
