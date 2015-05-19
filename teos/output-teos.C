// output-teos.x
// For a TEOS molecule simulation
// Quan (Andy) Wan Thu Nov 21 21:51:47 PST 2013

#include <iostream>
#include <iomanip>
#include <vector>
//#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "blas.h"
#include "D3vector.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "Stat.h"
#include "Tensor.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<6)
  {
    cout << "USAGE:\n";
    cout << "output-teos.x [input xyz] [input mlwf] [input quadrupole] [input polar] #frame\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad,fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fquad.open(argv[3]);
  fpolar.open(argv[4]);

  ofstream fdipole;
  fdipole.open("dipole");

  int nframe = atoi(argv[5]), nmo = 40, natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;



  for ( int iframe = 0; iframe < nframe; iframe ++)
  {


    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;
    //cout << "Frame" << iframe << endl;

    //read header
    tm_read.start();
    fxyz >> natom;
    fxyz >> c;
    tm_read.stop();  

    if ( iframe == 0 )
    {
      c.invert();
      atomset.resize(natom);
      mlwfset.resize(nmo);
    }

    //read atoms
    int icount = 0;
    tm_read.start();
    for ( int iatom = 0; iatom < natom; iatom ++ )
    {
      string name;

      fxyz >> name;

      if ( iframe == 0 )
      {
        if ( name.compare("O") == 0 ) 
        {
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6.0);
        }
        else if ( name.compare("Si") == 0 )
        {
          atomset[iatom].setmass(28.086);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(4.0);
        }
        else if ( name.compare("C") == 0 )
        {
          atomset[iatom].setmass(12.00);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(4.0);
        }
        else if ( name.compare("H") == 0 )
        {
          atomset[iatom].setmass(1.008);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(1.0);
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].x()<< endl;

    }//for iatom
    tm_read.stop();
    
    
    for ( int imo = 0; imo < nmo; imo++)
    {
      Mlwf * pwf = & mlwfset[imo];
      pwf -> readx(fmlwf);
      pwf -> reads(fmlwf);
      pwf -> readquad(fquad);
      pwf -> readp(fpolar);
      pwf -> setnumber(imo);
      //cout << imo << " " << pwf -> x() << endl;
    }

    D3vector dipole(0.0,0.0,0.0);

    for ( int i = 0; i < atomset.size(); i ++ )
      dipole += atomset[i].x() * atomset[i].charge();

    for ( int i = 0; i < mlwfset.size(); i ++ )
      dipole -= 2.0 * mlwfset[i].x();

    fdipole << dipole << endl;
  
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


  fxyz.close();
  fmlwf.close();
  fquad.close();

}//main
