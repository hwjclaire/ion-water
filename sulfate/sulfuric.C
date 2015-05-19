#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include "D3vector.h"
#include "Water.h"
#include "Sulfate.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#define AU2ANG 0.52918

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  tm.start();
  if (argc<3)
  {
    cout << "USAGE:\n";
    //cout << "sulfuric.x [input xyz] [input mlwf] #frame nMO nskip\n";
    cout << "sulfuric.x [input xyz] #frame\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf;
  fxyz.open(argv[1]);
  //fmlwf.open(argv[2]);

  int nhso4 = 0, nso4 = 0;


  ofstream fh3o, foh;
  fh3o.open("h3o.dat");
  foh.open("oh.dat");

  int nframe = atoi(argv[2]), natom;
  //int nmo = atoi(argv[4]), nskip = atoi(argv[5]);
  Cell c;

  for ( int iframe = 0; iframe < nframe; iframe ++)
  {

    std::vector<Atom> atomset;
    std::vector<Mlwf> mlwfset;
    std::vector<Water> waterset;


    //read header
    if ( iframe % 1000 == 0 ) 
      cout << "Frame" << iframe + 1 << endl;

    fxyz >> natom;
    fxyz >> c;

    if ( iframe == 0 ) c.invert();
    Sulfate sulfate(0,c);
    atomset.resize(natom);
    //mlwfset.resize(nmo);

    //read atoms
    int icount = 0;
    for ( int iatom = 0; iatom < natom; iatom ++ )
    { 
      string name;
      fxyz >> name;
      if ( 1 )
      {
        if ( name.compare("S") == 0 )
        {
          atomset[iatom].setmass(32.066);
          atomset[iatom].setname(name);
          sulfate.add_sulfur(atomset[iatom]);
        }
        if ( name.compare("K") == 0)
        {
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
          sulfate.add_oxygen(atomset[iatom]);
        }
        if ( name.compare("O") == 0 ) 
        {
          waterset.push_back( Water(icount ++ , c) );
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
        }
        if ( name.compare("H") == 0 || name.compare("Q") == 0 )
        {
          atomset[iatom].setmass(2.014);
          atomset[iatom].setname("H");
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].pos()<< endl;

    }//for iatom
  
    //assign H to O
    if ( 1 )  //assign molecules 
    {
      if ( iframe == 0)
      {
        cout << "natom " << atomset.size() << endl;
        cout << "nwater " << waterset.size() << endl;
      }

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
        Atom *patom = &atomset[iatom];
        if ( atomset[iatom].name().compare("H") ) continue;//is not H
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

      int noh3 = 0, noh = 0;
      double soh3[3], soh = 0;
      soh3[0] = soh3[1] = soh3[2] = 0;
      sulfate.check_atom();
      
      for ( vector<Water>::iterator i = waterset.begin(); i < waterset.end(); i ++ )
      {
        //i -> check_atom();
        if ( i -> isoh() )
        {
          noh ++;
          soh = sulfate.distance(&(*i));
          //cout << "oh" << endl;
        }
        if ( i -> isoh3() )
        {
          soh3[noh3++] = sulfate.distance(&(*i));

        }//if isoh3

      }//for i

      if( noh3 == 2 && soh3[0] > soh3[1] ) 
      {
        double tmp;
        tmp = soh3[0];
        soh3[0] = soh3[1] ;
        soh3[1] = tmp;
      }
      else if ( noh3 == 3 )
      {
        double mid = soh3[0], max = mid, min = mid;

        if ( soh3[2] > max ) max = soh3[2];
        if ( soh3[2] < min ) min = soh3[2];
       
        if ( soh3[1] > max ) max = soh3[1];
        if ( soh3[1] < min ) min = soh3[1];
    
        for ( int i = 0; i < 3; i ++)
        {
          if ( soh3[i] > min && soh3[i] < max )
          {
            mid = soh3[i];
            break;
          }
        }
        soh3[0] = min; soh3[1] = mid; soh3[2] = max;
      }

      fh3o.width(10);
      foh.width(10);

      fh3o  << iframe +1 << "  " << noh3 << " " << noh << " ";
      for ( int i = 0; i < 3 - noh3 ; i ++)  fh3o << setw(10) << 0.0;
      for ( int i = 0; i < noh3 ; i ++) fh3o << "  " << setw(10) <<  soh3[i] * AU2ANG ;
      fh3o << endl;

      if ( noh == 1 ) foh << iframe + 1 << "  "  << soh * AU2ANG << endl;

          
    }// if 1

    if ( sulfate.minus1() ) nhso4 ++;
    if ( sulfate.minus2() ) nso4 ++;

    
     
#if 0
    if ( iframe % nskip == 0 )
    {
      //add MLWFs to water
      for ( int i = 0; i < waterset.size(); i ++ ) waterset[i].reset_wf();
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        //pwf -> readx(fmlwf);
        //pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x << endl;

        double min_dist = 100;
        Mol * min;
        int nwater = waterset.size();
        for ( int iwater = 0; iwater < nwater; iwater ++)
        {
          Water * pwater = &waterset[iwater];
          //if (waterset[iwater].wf_full()) continue;
          double dist = pwater -> distance(pwf,c);
          if ( dist < min_dist )
          {
            min = pwater;
            min_dist = dist;
          }
        }
        min -> add_wf(*pwf, c);
        //cout << imo << " " << min_dist << endl;
      }
  
      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].check_wf();
      }

      //output
  
    }// if iframe % nskip
#endif


  




  }// for iframe
  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  cout << nhso4 << "   " << nso4 << endl;


  fxyz.close();
  //fmlwf.close();
  fh3o.close();
  foh.close();

}//main

