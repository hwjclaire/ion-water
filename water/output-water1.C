// output-spread.x
// Given a d2o simulation including polarizability
// output coordinates of individual waters in qbox input format
// Quan (Andy) Wan Thu Dec 13 14:23:19 PST 2012

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
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

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<9)
  {
    cout << "USAGE:\n";
    cout << "output-water1.x [input xyz] [NULL] [NULL]  #frame nmo nskip1 [NULL] nwater\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fpolar;
  fxyz.open(argv[1]);
  //fmlwf.open(argv[2]);
  //fpolar.open(argv[3]);



  int nframe = atoi(argv[4]), nmo = atoi(argv[5]), nskip1 = atoi(argv[6]), natom;
  int nwater = atoi(argv[8]), nskip2 = atoi(argv[7]);
  assert ( nskip1 > 0 );
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  ofstream polarmlwf[nwater];
  ofstream head;
  head.open("head");


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
        if ( name.compare("H") == 0 )
        {
          atomset[iatom].setmass(2.014);
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
        name << "waterCoord" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmlwf[i].open(fptr);
        polarmlwf[i].setf(ios::fixed, ios::floatfield);
        polarmlwf[i].setf(ios::right, ios::adjustfield);
        polarmlwf[i].precision(4);
      }

    }// if ifram == 0;  
    //cout << "addh\n";

    //if ( nskip1 == 0 || iframe % nskip1 == 0 )
    if (0)
    {

      tm_read.start();     
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        //pwf -> readp(fpolar);
        //pwf -> setnumber(imo);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();


      if ( iframe % ( nskip1 * nskip2 ) != 0 ) continue;
  

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
      {
        //waterset[i].mlwfmove();
        waterset[i].check_wf();
        //waterset[i].Compute_cc();
        //waterset[i].Compute_cm();
        //waterset[i].Compute_dipole();
      }
      tm_misc.stop();
      //cout << "Check Done" << endl;
  

      tm_output.start();


      tm_output.stop();

    }// if 0


    if ( iframe % nskip1 == 0 ) 
      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        if ( iframe == 0 )
        {
          if ( iwater == 0 )
          {
            head  << "atom O1 oxygen" << waterset[iwater].o()->x() << endl
                  << "atom H1 hydrogen" << waterset[iwater].h(0)->x() << endl
                  << "atom H2 hydrogen" << waterset[iwater].h(1)->x() << endl
                  << "set ecut 20\nrun 0 50\n"
                  << "set ecut 40\nrun 0 30\n"
                  << "set ecut 60\nrun 0 30\n"
                  << "set ecut 75\nrun 0 20\n"
                  << "set ecut 85\nrun 0 20\n"
                  << "save water" << iwater+1 << ".xml\n";
          }
          else
          {
            head  << "move O1 to" << waterset[iwater].o()->x() << endl
                  << "move H1 to" << waterset[iwater].h(0)->x() << endl
                  << "move H2 to" << waterset[iwater].h(1)->x() << endl
                  << "set ecut 0\n"
                  << "set ecut 20\nrun 0 50\n"
                  << "set ecut 40\nrun 0 30\n"
                  << "set ecut 60\nrun 0 30\n"
                  << "set ecut 75\nrun 0 20\n"
                  << "set ecut 85\nrun 0 20\n"
                  << "save water" << iwater+1 << ".xml\n";
          }                      
        }
        else
        {
          polarmlwf [ iwater ] << "move O1 to" << waterset[iwater].o()->x() << endl;
          polarmlwf [ iwater ] << "move H1 to" << waterset[iwater].h(0)->x() << endl;
          polarmlwf [ iwater ] << "move H2 to" << waterset[iwater].h(1)->x() << endl;
          polarmlwf [ iwater ] << "run 0 30\n";
        }

      }
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
  fpolar.close();

}//main
