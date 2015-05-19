// conver Cui's pos files to au XYZ files
// Given a d2o simulation 
// Quan (Andy) Wan Thu Dec 13 14:23:19 PST 2012

#include <iostream>
#include <iomanip>
#include <vector>
//#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "Timer.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<6)
  {
    cout << "USAGE:\n";
    cout << "cui2andy.x [input O] [ input H ] [output] #frame nwater\n";
    return 1;
  }
 
  ifstream fxyzo, fxyzh;
  fxyzo.open(argv[1]);
  fxyzh.open(argv[2]);
  
  ofstream fxyz;
  fxyz.open(argv[3]);
  fxyz.precision(9);


  string str_O = "O", str_H = "H";

  int nframe = atoi(argv[4]);
  int nwater = atoi(argv[5]), natom = nwater * 3;
  double volume;

  std::vector < double > cell;
  cell.resize(9,0.0);
  if ( nwater == 64 )
  {
    cell[0] = 23.464;
    cell[4] = 23.464;
    cell[8] = 23.464;
  }
  
  if ( nwater == 32 )
  {
    cell[0] = 18.637050;
    cell[4] = 18.637050;
    cell[8] = 18.637050;
  }


  for ( int iframe = 0; iframe < nframe; iframe ++)
  {

    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;

    fxyz << natom << endl;
    for ( int i = 0; i < 9; i ++ ) fxyz << cell[i]<< " ";
    fxyz << endl;
  

    //read atoms
    int icount = 0;
    tm_read.start();
    for ( int iatom = 0; iatom < nwater*3; iatom ++ )
    {
      double t[3];
      if ( iatom < nwater ) 
      {
        for ( int i = 0; i < 3; i ++ ) fxyzo >> t[i];
        fxyz << str_O << " ";
        for ( int i = 0; i < 3; i ++ ) fxyz << t[i] << " ";
        for ( int i = 0; i < 3; i ++ ) fxyz << 0.0 << " ";
        fxyz << endl;
      }
      else
      {
        for ( int i = 0; i < 3; i ++ ) fxyzh >> t[i];
        fxyz << str_H << " ";
        for ( int i = 0; i < 3; i ++ ) fxyz << t[i] << " ";
        for ( int i = 0; i < 3; i ++ ) fxyz << 0.0 << " ";
        fxyz << endl;
      }

    }//for iatom
    tm_read.stop();


  }// for iframe

  tm.stop();
  cout .precision(3);
  cout << "ALL:    Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;
  cout << "Read:   Real Time: " << tm_read.real() <<  "s  CPU Time:" << tm_read.cpu() << "s" << endl;
  cout << "Addwf:  Real Time: " << tm_addwf.real() <<  "s  CPU Time:" << tm_addwf.cpu() << "s" << endl;
  cout << "Alpha:  Real Time: " << tm_alpha.real() <<  "s  CPU Time:" << tm_alpha.cpu() << "s" << endl;
  cout << "Rot:    Real Time: " << tm_rot.real() <<  "s  CPU Time:" << tm_rot.cpu() << "s" << endl;
  cout << "output: Real Time: " << tm_output.real() <<  "s  CPU Time:" << tm_output.cpu() << "s" << endl;
  cout << "Misc:   Real Time: " << tm_misc.real() <<  "s  CPU Time:" << tm_misc.cpu() << "s" << endl;


  fxyzo.close();
  fxyzh.close();
  fxyz.close();

}//main
