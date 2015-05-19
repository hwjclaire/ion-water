// output-spread-forCui.x //for Cui's data files 
// Given a d2o simulation including polarizability
// output sum of MLWF spread all waters
// Quan (Andy) Wan Thu Dec 13 14:23:19 PST 2012

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

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<9)
  {
    cout << "USAGE:\n";
    cout << "output-spread-forCui.x [input O] [ input H ] [input mlwf] [ input spread ]  #frame nMO nskip nwater\n";
    return 1;
  }
 
  ifstream fxyzo, fxyzh, fmlwf, fspread;
  fxyzo.open(argv[1]);
  fxyzh.open(argv[2]);
  
  fmlwf.open(argv[3]);
  fspread.open(argv[4]);

  string str_O = "O", str_H = "H";

  int nframe = atoi(argv[5]), nmo = atoi(argv[6]), nskip = atoi(argv[7]);
  int nwater = atoi(argv[8]), natom = nwater * 3;
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

  Cell c(cell);

  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  ofstream polarmlwf[nwater];


  for ( int iframe = 0; iframe < nframe; iframe ++)
  {


    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;

    //read header

    if ( iframe == 0 )
    {
      c.invert();
      atomset.resize(natom);
      mlwfset.resize(nmo);
    }

    //read atoms
    int icount = 0;
    tm_read.start();
    for ( int iatom = 0; iatom < nwater*3; iatom ++ )
    {
      //cout << iatom << endl;
      if ( iframe == 0 )
      {
        if ( iatom < nwater ) 
        {
          waterset.push_back( Water(icount ++, c ) );
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(str_O);
          atomset[iatom].setcharge(6.0);
        }
        else
        {
          atomset[iatom].setmass(2.014);
          atomset[iatom].setname(str_H);
          atomset[iatom].setcharge(1.0);
        }
      }//if iframe == 0
      // read position, velocity
      if ( iatom < nwater ) 
        atomset[iatom].readx(fxyzo);
      else
        atomset[iatom].readx(fxyzh);
      //cout << atomset[iatom].x()<< endl;

    }//for iatom
    tm_read.stop();

    //cout << "Read\n";
  
    //assign H to O
    if ( iframe == 0 )  //assign molecules 
    {

      cout << "natom " << atomset.size() << endl;
      cout << "nwater " << waterset.size() << endl;

      int watercount = 0;
      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        if ( atomset[iatom].name().compare("O")  ) continue;//is not O
        //cout << watercount << " " << & waterset[watercount] << " " << & atomset[iatom] << endl;
        waterset[watercount].add_oxygen ( atomset[iatom] );
        watercount ++;
      }
      //assign hydrogen to oxygen
      for ( int iatom = 0; iatom < natom; iatom ++)
      { 
        //cout << iatom << endl;
        Atom * patom = & atomset[iatom];
        if ( patom -> name().compare("H") ) continue;//is not H
        double min_dist = 1000; Water * min;
        for ( int iwater = 0; iwater < waterset.size(); iwater ++)
        {
          //if ( pwater -> atom_full ) continue; //water full
          Water * pwater = &waterset[iwater];
          double dist = pwater -> distance(patom);
  
          if ( dist < min_dist) 
          {
            min = pwater;
            min_dist = dist;
          }//if

        }//for jatom
        min -> add_hydrogen ( *patom ); 
      }

      //cout << "addH\n";

      assert ( waterset.size() * 4 == nmo );
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        waterset[i].check_atom();
      }

      //cout << "check\n";
      
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        ostringstream name;
        name << "waterspread" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmlwf[i].open(fptr);
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
        pwf -> readpos(fmlwf);
        pwf -> reads(fspread);
        //pwf -> readp(fpolar);
        //pwf -> setnumber(imo);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();


      if ( nskip != 0 && iframe % ( nskip * 4 ) != 0 ) continue;
  

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

#if 1
      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        double sum = 0;
        for ( int imlwf = 0; imlwf < 4; imlwf ++)
          polarmlwf [ iwater ] << setw(12) << 
                waterset[iwater].wf(imlwf)->spread();


        polarmlwf [ iwater ] << setw(12) << 
                waterset[iwater].wf(0)->spread()
              + waterset[iwater].wf(1)->spread();

        polarmlwf [ iwater ] << setw(12) <<  
                waterset[iwater].wf(2)->spread()
              + waterset[iwater].wf(3)->spread();

        polarmlwf [ iwater ] << setw(12) <<
                waterset[iwater].wf(0)->spread()
              + waterset[iwater].wf(1)->spread()
              + waterset[iwater].wf(2)->spread()
              + waterset[iwater].wf(3)->spread();


        polarmlwf [ iwater ] << endl;

      }//for iwater

#endif

      tm_output.stop();

    }// if iframe % nskip

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
  fmlwf.close();
  fspread.close();

}//main
