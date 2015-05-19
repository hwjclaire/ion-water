#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "D3vector.h"
#include "Water.h"
#include "Sulfate.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "Stat.h"
#define AU2ANG 0.52918
#define rshell 4.45

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  tm.start();
  if (argc<7)
  {
    cout << "USAGE:\n";
    //cout << "sulfuric.x [input xyz] [input mlwf] #frame nMO nskip\n";
    cout << "sulfate-spread.x [input xyz] [input mlwf] #frame nMO nskip nskip2\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);

  

  int nframe = atoi(argv[3]), natom;
  int nmo = atoi(argv[4]), nskip = atoi(argv[5]), nskip2 = atoi(argv[6]);
  Cell c;

  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  mlwfset.resize(nmo);
  std::vector<Water> waterset;
  Sulfate sulfate(0,c);


  ofstream fwater[62],fsulfate;
  for ( int i = 0; i < 62; i ++ )
  {
    ostringstream name;
    name << "waterspread" << i + 1 << ".dat";
    string filename = name.str();
    char * fptr = (char*)filename.c_str();
    fwater[i].open(fptr);
    fwater[i].setf(ios::fixed, ios::floatfield);
    fwater[i].setf(ios::right, ios::adjustfield);
    fwater[i].precision(4);
  }

  fsulfate.open("sulfatespread.dat");
  fsulfate.setf(ios::fixed, ios::floatfield);
  fsulfate.setf(ios::right, ios::adjustfield);
  fsulfate.precision(4);
  

  for ( int iframe = 0; iframe < nframe; iframe ++)
  {

    //read header
    if ( iframe % 1000 == 0) cout << "Frame" << iframe  << endl;
    //cout << "Frame" << iframe  << endl;
    fxyz >> natom;
    fxyz >> c;

    if ( iframe == 0 ) c.invert();
    atomset.resize(natom);
    //mlwfset.resize(nmo);

    //read atoms
    int icount = 0;
    for ( int iatom = 0; iatom < natom; iatom ++ )
    { 
      string name;
      fxyz >> name;
      if ( iframe == 0 )
      {
        if ( name.compare("S") == 0 )
        {
          atomset[iatom].setmass(32.066);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6);
          sulfate.add_sulfur(atomset[iatom]);
        }
        if ( name.compare("K") == 0)
        {
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6);
          sulfate.add_oxygen(atomset[iatom]);
        }
        if ( name.compare("O") == 0 ) 
        {
          waterset.push_back( Water(icount ++,c ) );
          atomset[iatom].setmass(15.999);
          atomset[iatom].setcharge(6);
          atomset[iatom].setname(name);
        }
        if ( name.compare("H") == 0 || name.compare("Q") == 0 )
        {
          atomset[iatom].setmass(2.014);
          atomset[iatom].setcharge(1);
          atomset[iatom].setname("H");
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].pos()<< endl;

    }//for iatom

    //cout << "ReadAtom\n";
  
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
        if ( atomset[iatom].name().compare("H") ) continue;//is not H
        Atom *patom = &atomset[iatom];
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
        double distsh = sulfate.distance(patom);
        if ( distsh < min_dist)
        {
          sulfate.add_hydrogen( *patom );
        }
        else
        {
          min -> add_hydrogen ( *patom ); 
        }
      }

      sulfate.check_atom();
      
    }// if ifram == 0;  
     
    //cout << "Addh\n";

    //MLWF
    if ( nskip == 0 || iframe % nskip == 0 )
    {
      //add MLWFs to water
      if ( nskip2 == 0 || iframe %( nskip * nskip2) == 0 )
      {
        sulfate.reset_wf();
        for ( int i = 0; i < waterset.size(); i ++ )
        {
          waterset[i].watermove();
          waterset[i].reset_wf(); 
        }
      }


      for ( int imo = 0; imo < nmo; imo++)
      {
        //Read from file
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        //pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x() << "  ";

        if ( nskip2 != 0 && iframe %( nskip * nskip2) != 0 ) continue;
        
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
        //cout << imo << " " << min_dist << "  ";
        double sh = sulfate.distance(pwf);
        if ( min_dist < sh )
        {
          min -> add_wf(*pwf);
        }
        else
        {
          sulfate.add_wf(*pwf);
        }
      }//for imo 

      //cout << "Addmlwf\n";
      if ( nskip2 != 0 && iframe %( nskip * nskip2) != 0 ) continue;


      double sspread = 0; // sum of sulfate spread;
      for ( int i = 0; i < 16; i++)
      {
        sspread += sulfate.wf(i)->spread();
      }

      fsulfate << sspread << endl;


      double wspreadtot = 0;
      for ( int i = 0; i < 62; i ++ )
      {
        double wspread = 0; //sum of water spread 
        for ( int ii = 0; ii < 4; ii++)  
        {
          wspread += waterset[i] . wf(ii) -> spread();
          fwater[i] << setw(12) << waterset[i] . wf(ii) -> spread();
        }//for ii
        wspreadtot += wspread;

        fwater[i] << setw(12) << waterset[i] . wf(0) -> spread()
                               + waterset[i] . wf(1) -> spread();

        fwater[i] << setw(12) << waterset[i] . wf(2) -> spread()
                               + waterset[i] . wf(3) -> spread();

        fwater[i] << setw(12) << wspread << endl;

      }//for i
      

    }// if iframe % nskip


  }// for iframe
  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  fxyz.close();
  fsulfate.close();
  
  for ( int i = 0; i < 62; i ++ )
    fwater[i].close();

}//main

