// pick-inter.-fs2x
// Given a d2o simulation including polarizability
// pick the molecules that are in 2nd solvation shell of water
//
// Example:
// [step]
//  [water1]
//    [number of water in the second shell of water1] [water1] [water2] [water3] ....
//  [water2]
//    [number of water in the second shell of water1] [water1] [water2] [water3] ....
//  [water3]
//    [number of water in the second shell of water1] [water1] [water2] [water3] ....
//  ...
// [step]
//  ...
// ...
// 
// Caveat: water count 1-64
//
// Andy Wan  Mon Dec 10 03:04:44 EST 2012

#include <iostream>
#include <iomanip>
#include <vector>
//#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "D3vector.h"
#include "Mol.h"
#include "Water.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "Stat.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<8)
  {
    cout << "USAGE:\n";
    cout << "pickNB.x [input xyz] [input mlwf] [input polariz]  #frame nMO nskip [output]\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fpolar.open(argv[3]);

  ofstream hbond,fhbstat;
  hbond.open(argv[7]);  // output stat
  fhbstat.open("hbond.dat");  // output stat about second shell

  int nframe = atoi(argv[4]), nmo = atoi(argv[5]), nskip = atoi(argv[6]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  Stat hbstat(0,50);

  for ( int iframe = 0; iframe < nframe; iframe ++)
  {

    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;

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
          waterset.push_back( Water(icount ++, c ) );
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6.0);
        }
        else
        {
          atomset[iatom].setmass(2.014);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(1.0);
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].pos()<< endl;

    }//for iatom
    tm_read.stop();

    //assign H to O
    if ( iframe == 0 )  //assign molecules 
    {

      cout << "natom " << atomset.size() << endl;
      cout << "nwater " << waterset.size() << endl;

      int watercount = 0;
      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        if ( atomset[iatom].name().compare("O")  ) continue;//is not O
        //cout << watercount << & waterset[watercount] << " " << & atomset[iatom] << endl;
        waterset[watercount].add_oxygen ( atomset[iatom] );
        watercount ++;
      }
      //assign hydrogen to oxygen
      for ( int iatom = 0; iatom < natom; iatom ++)
      { 
        Atom * patom = & atomset[iatom];
        if ( patom -> name().compare("H") ) continue;//is not H
        double min_dist = 1000; Mol * min;
        for ( int iwater = 0; iwater < waterset.size(); iwater ++)
        {
          //if ( pwater -> atom_full ) continue; //water full
          Mol* pwater = &waterset[iwater];
          double dist = pwater -> distance(patom);

          if ( dist < min_dist) 
          {
            min = pwater;
            min_dist = dist;
          }//if

        }//for jatom
        min -> add_hydrogen ( *patom ); 
      }

      for ( int i = 0; i < waterset.size(); i ++ )
      {
        waterset[i].check_atom();
      }

    }// if ifram == 0;  
    //cout << "addh\n";

    if ( nskip == 0 || iframe % nskip == 0 )
    {

      tm_read.start();     
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();

      if ( iframe % ( nskip * 4 ) != 0 ) continue;
  

      tm_misc.start();
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        //update oh distances.
        waterset[i].watermove();

        waterset[i].reset_wf();
      }
      tm_misc.stop();
      //cout << "resetwf\n";

      tm_addwf.start();
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        double min_dist = 100;
        Mol * min;
        int nwater = waterset.size();
        for ( int iwater = 0; iwater < nwater; iwater ++)
        {
          Water * pwater = &waterset[iwater];
          //if (waterset[iwater].wf_full()) continue;
          double dist = pwater -> distance(pwf);
          if ( dist < min_dist )
          {
            min = pwater;
            min_dist = dist;
          }
        }
        min -> add_wf(*pwf);
        //cout << imo << " " << min_dist << endl;
      }
      tm_addwf.start();
      //cout << "addwf\n"; 


      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].check_wf();
      }
      //cout << "Check Done" << endl;
  

      hbond << waterset.size() << endl;
      
      const double thresh1 = 3.35 / 0.52918;
      const double thresh2 = 5.55 / 0.52918;
      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {

        tm_misc.start();

        std::vector < int > secondshell;

        const D3vector& vec1 = waterset[iwater].o() -> x();//o position

        for ( int jwater = 0; jwater < waterset.size(); jwater ++ )
        {

          if ( iwater == jwater ) continue;
          D3vector oo = waterset[jwater].o() -> x() - vec1;
          c.images(oo);
          double d_oo = length(oo);

          //cout << iwater << "    " << jwater  << "   " << oo << "  " << d_oo << "  ";

          if ( d_oo > thresh1 && d_oo < thresh2)
            secondshell.push_back(jwater+1);            

          //cout << endl;
        }//for jwater

        tm_misc.stop();

        tm_output.start();
      
        //output second solvation shell
        hbond << "    " << secondshell.size() << "  ";
        for ( int ifs = 0; ifs < secondshell.size(); ifs ++ )
          hbond << secondshell[ifs] << "   ";
        hbond << endl;

        tm_output.stop();
  
        hbstat.add(secondshell.size());

      }//for iwater

    }// if iframe % nskip

  }// for iframe

  hbstat.print(fhbstat);


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
  fpolar.close();

}//main
