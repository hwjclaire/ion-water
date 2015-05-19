// pickNB-lp.x
// Given a d2o simulation including polarizability
// pick the molecules that at least one of its OH bond is not hydrogen bonded
// Andy Wan Mon Dec 10 03:05:41 EST 2012

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
  if (argc<7)
  {
    cout << "USAGE:\n";
    cout << "pickNB-lp.x [input xyz] [input mlwf] [input polariz]  #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fpolar.open(argv[3]);


  ofstream polarmlwfaxis[128];


  int nframe = atoi(argv[4]), nmo = atoi(argv[5]), nskip = atoi(argv[6]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

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

      for ( int i = 0; i < 128; i ++ )
      {
        ostringstream name;
        name << "polarmlwf" << i+1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmlwfaxis[i].open(fptr);
        polarmlwfaxis[i].setf(ios::fixed, ios::floatfield);
        polarmlwfaxis[i].setf(ios::right, ios::adjustfield);
        polarmlwfaxis[i].precision(4);
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
        pwf -> setnumber(imo);
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


      tm_misc.start();

      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].mlwfmove();
        waterset[i].check_wf();
        //cout << i << "  " ;
        //for ( int j = 0; j < 4; j ++ ) cout << waterset[i].wf(j) -> number() << "  ";
        //cout << endl;
      }
      tm_misc.stop();
      //cout << "Check Done" << endl;
  

      tm_output.start();


      double thresh1 = 3.35 / 0.52918;
      double thresh2 = sqrt(3) / 2;
      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        const D3vector& vec1 = waterset[iwater].o() -> x();//o position

        D3vector obp[2];//O-LP vector 
        obp[0]  = normalized ( waterset[iwater].owf(2) );
        obp[1]  = normalized ( waterset[iwater].owf(3) );

        bool hb[2]; hb[0] = false; hb[1] = false;
        //a oh is hbonded

        for ( int jwater = 0; jwater < waterset.size(); jwater ++ )
        {

          if ( iwater == jwater ) continue;
          if ( hb[0] == true && hb[1] == true ) continue;
          D3vector oo = waterset[jwater].o() -> x() - vec1;
          c.images(oo);
          double d_oo = length(oo);

          //cout << iwater << "    " << jwater  << "   " << oo << "  " << d_oo << "  ";

          if ( d_oo < thresh1 )
          {
            oo /= d_oo; // oo is normalized;

            for ( int ioh = 0; ioh < 2; ioh ++ )
            {
              if ( obp[ioh] * oo > thresh2 ) hb[ioh] = true;
              //cout << oh[ioh] * oo << "   ";
            }
          }
          //cout << endl;
        }//for jwater

#if 1
        for ( int ioh = 0; ioh < 2; ioh ++ )
        {
          if ( hb[ioh] )
            polarmlwfaxis [ iwater * 2 + ioh ] << "1";
          else
            polarmlwfaxis [ iwater * 2 + ioh ] << "0";
        }

        polarmlwfaxis [ iwater * 2 ] 
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[0]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[4]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[8]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[1]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[2]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[5]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[3]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[6]
                          << setw(8) << waterset[iwater].wf(2)->polar_tensor()[7]
                          << endl;

        polarmlwfaxis [ iwater * 2 + 1 ] 
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[0]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[4]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[8]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[1]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[2]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[5]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[3]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[6]
                          << setw(8) << waterset[iwater].wf(3)->polar_tensor()[7]
                          << endl;
#else
        if ( hb[0] and hb[1] )
        {
          polarmlwfaxis [ iwater * 2 ] << "1";
          polarmlwfaxis [ iwater * 2 + 1 ] << "1";
        }
        else
        {
          polarmlwfaxis [ iwater * 2 ] << "0";
          polarmlwfaxis [ iwater * 2 + 1 ] << "0";
        }

        polarmlwfaxis [ iwater * 2 ]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[0]
                                      + waterset[iwater].wf(1)->polar_tensor()[0]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[4]
                                      + waterset[iwater].wf(1)->polar_tensor()[4]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[8]
                                      + waterset[iwater].wf(1)->polar_tensor()[8]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[1]
                                      + waterset[iwater].wf(1)->polar_tensor()[1]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[2]
                                      + waterset[iwater].wf(1)->polar_tensor()[2]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[5]
                                      + waterset[iwater].wf(1)->polar_tensor()[5]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[3]
                                      + waterset[iwater].wf(1)->polar_tensor()[3]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[6]
                                      + waterset[iwater].wf(1)->polar_tensor()[6]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[7]
                                      + waterset[iwater].wf(1)->polar_tensor()[7]
                          << endl;

        polarmlwfaxis [ iwater * 2 + 1 ]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[0]
                                      + waterset[iwater].wf(1)->polar_tensor()[0]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[4]
                                      + waterset[iwater].wf(1)->polar_tensor()[4]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[8]
                                      + waterset[iwater].wf(1)->polar_tensor()[8]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[1]
                                      + waterset[iwater].wf(1)->polar_tensor()[1]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[2]
                                      + waterset[iwater].wf(1)->polar_tensor()[2]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[5]
                                      + waterset[iwater].wf(1)->polar_tensor()[5]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[3]
                                      + waterset[iwater].wf(1)->polar_tensor()[3]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[6]
                                      + waterset[iwater].wf(1)->polar_tensor()[6]
                          << setw(8) << waterset[iwater].wf(0)->polar_tensor()[7]
                                      + waterset[iwater].wf(1)->polar_tensor()[7]
                          << endl;

#endif

      }//for iwater

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


  fxyz.close();
  fmlwf.close();
  fpolar.close();

}//main
