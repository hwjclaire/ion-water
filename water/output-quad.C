// output-quad.x
// Given a d2o simulation computing quadrupole of MLWFs, output quadrupoles
// Quan (Andy) Wan Fri Sep 13 15:20:33 PDT 2013

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
#include "Tensor.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<7)
  {
    cout << "USAGE:\n";
    cout << "output-quad.x [input xyz] [input mlwf] [input quadrupole]  #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fquad.open(argv[3]);



  int nframe = atoi(argv[4]), nmo = atoi(argv[5]), nskip = atoi(argv[6]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  ofstream polarmlwf[64];

  Stat quad1(-4,4,0.02);
  Stat quad2(-4,4,0.02);
  Stat quad3(-4,4,0.02);
  Stat dip(0,2,0.02);

  ofstream fdip,fq1,fq2,fq3;
  fdip.open("statdip");
  fq1.open("statquad1");
  fq2.open("statquad2");
  fq3.open("statquad3");

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
          atomset[iatom].setmass(1.0079);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(1.0);
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].x()<< endl;

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


      assert ( waterset.size() * 4 == nmo );
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        waterset[i].check_atom();
      }

      
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        ostringstream name;
        name << "waterquad" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        //polarmlwf[i].open(fptr);
        polarmlwf[i].setf(ios::fixed, ios::floatfield);
        polarmlwf[i].setf(ios::right, ios::adjustfield);
        polarmlwf[i].precision(4);
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
        pwf -> reads(fmlwf);
        pwf -> readquad(fquad);
        //pwf -> setnumber(imo);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();


      if ( iframe % ( nskip ) != 0 ) continue;
  

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

      tm_misc.start();
      for ( int i = 0; i < waterset.size(); i ++ )
        waterset[i].mlwfmove();

      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].check_wf();
        waterset[i].Compute_cc();
        waterset[i].Compute_cm();
        waterset[i].Compute_dipole();
        waterset[i].Compute_quad();
        quad1.add ( waterset[i].quad_eigval()[0] );
        quad2.add ( waterset[i].quad_eigval()[1] );
        quad3.add ( waterset[i].quad_eigval()[2] );
        dip.add(length(waterset[i].dipole()));
      }
      tm_misc.stop();
      //cout << "Check Done" << endl;
  
      tm_output.start();

#if 1
      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
        polarmlwf [ iwater ] << waterset[iwater].quad() << endl;

#endif
      tm_output.stop();

    }// if iframe % nskip

  }// for iframe

  quad1.print(fq1);
  quad2.print(fq2);
  quad3.print(fq3);
  dip.print(fdip);
  cout << dip.avg() << endl;
  cout << quad1.avg() << endl;
  cout << quad2.avg() << endl;
  cout << quad3.avg() << endl;

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
