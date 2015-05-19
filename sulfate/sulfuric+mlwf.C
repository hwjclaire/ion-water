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
#include "Stat.h"
#define AU2ANG 0.52918

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  tm.start();
  if (argc<6)
  {
    cout << "USAGE:\n";
    cout << "sulfuric+mlwf.x [input xyz] [input mlwf] #frame nMO nskip\n";
    //cout << "sulfuric.x [input xyz] #frame\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);

  ofstream fsulfatemlwfspread, fwatermlwfspread, fsulfatedipole, fwaterdipole;
  fsulfatemlwfspread.open("smlwfspread.dat");
  fwatermlwfspread.open("wmlwfspread.dat");
  fsulfatedipole.open("sdipole.dat");
  fwaterdipole.open("wdipole.dat");


  ofstream fsulfatemlwfdist, fwatermlwfdist, fsulfatesodist;
  fsulfatemlwfdist.open("smlwfdist.dat");
  fwatermlwfdist.open("wmlwfdist.dat");
  fsulfatesodist.open("sodist.dat");


  Stat sulfatemlwfspread(1,1.8);
  Stat watermlwfspread(1,1.8);
  Stat sulfatedipole(0,4.5);
  Stat waterdipole(0,10);
  Stat sulfatemlwfdist(1,4);
  Stat watermlwfdist(0.0,1.6);
  Stat sulfatesodist(2.2,4.5); 

  int nhso4 = 0, nso4 = 0;



  int nframe = atoi(argv[3]), natom;
  int nmo = atoi(argv[4]), nskip = atoi(argv[5]);
  Cell c;

  for ( int iframe = 0; iframe < nframe; iframe ++)
  {

    std::vector<Atom> atomset;
    std::vector<Mlwf> mlwfset;
    std::vector<Water> waterset;


    //read header
    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;
    fxyz >> natom;
    fxyz >> c;

    if ( iframe == 0 ) c.invert();
    Sulfate sulfate(0,c);
    atomset.resize(natom);
    mlwfset.resize(nmo);

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
          waterset.push_back( Water(icount ++ , c) );
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

    }// if 1

    sulfate.check_atom();
    for ( vector<Water>::iterator i = waterset.begin(); i < waterset.end(); i ++ )
      i -> check_atom();

    if ( sulfate.minus1() ) nhso4 ++;
    if ( sulfate.minus2() ) nso4 ++;  
     
    if ( iframe % nskip == 0 )
    {
      //add MLWFs to water
      sulfate.reset_wf();
      for ( int i = 0; i < waterset.size(); i ++ ) waterset[i].reset_wf();
      for ( int imo = 0; imo < nmo; imo++)
      {
        //Read from file
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        //pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x() << "  ";

        double min_dist = 100;
        Mol * min;
        int nwater = waterset.size();
        for ( std::vector <Water> :: iterator iwater = waterset.begin(); iwater < waterset.end(); iwater ++)
        {
          //if (waterset[iwater].wf_full()) continue;
          double dist = iwater -> distance(pwf);
          if ( dist < min_dist )
          {
            min = &(*iwater);
            min_dist = dist;
          }
        }
        //cout << imo << " " << min_dist << "  ";
        double sh = sulfate.distance(pwf);
        if ( min_dist < sh )
        {
          min -> add_wf(*pwf,true);
        }
        else
        {
          sulfate.add_wf(*pwf,true);
        }
      }//for imo 
      
      //cout << sulfate.nwf() <<endl;
      sulfate.check_wf();
      sulfate.Compute_cm();
      sulfate.Compute_dipole();
      //cout << sulfate.tot_charge() << "  " << sulfate.cm() << "  " << sulfate.dipole() << endl;

      //cout << sulfate.dipole() << endl;
      sulfatedipole.add ( length ( sulfate.dipole() ) );

      sulfatesodist.reset();
      for ( int i = 0; i < 4; i++ )
      {
        sulfatesodist.add ( sulfate.s()->distance(sulfate.o(i), &c) );
      }
  
      fsulfatesodist << iframe << "  " << sulfatesodist.min() 
                     << "  " << sulfatesodist.max() << "  "
                     << sulfatesodist.max() - sulfatesodist.min() <<  endl;

      for ( int i = 0; i < 16; i++)
      {
        sulfatemlwfspread.add(sulfate.wf(i)->spread());
        sulfatemlwfdist.add(sulfate.s()->distance(sulfate.wf(i), &c)); 
      }

      

      for ( vector<Water>::iterator i = waterset.begin(); i < waterset.end(); i ++ )
      {
        //i -> check_wf();
        i -> Compute_cm();
        i -> Compute_dipole();

        //cout << i -> o() -> x() << " " << i -> h(0) -> x() << " " << i -> h(1) -> x() << endl;
        //cout << i->tot_charge() << "  " << i->cm() << "  " << i->dipole() << endl;
        //cout << i -> dipole() << endl;
        waterdipole.add ( length( i -> dipole() ) );
        for ( int ii = 0; ii < 4; ii++)  
        {
          watermlwfspread.add ( i-> wf(ii) -> spread() );
          watermlwfdist.add ( i-> o() -> distance( i -> wf(ii), &c ) );
        }//for ii

      }//for i

    }// if iframe % nskip

  }// for iframe

  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  sulfatemlwfspread.print(fsulfatemlwfspread);
  watermlwfspread.print(fwatermlwfspread);
  sulfatedipole.print(fsulfatedipole);
  waterdipole.print(fwaterdipole);
  sulfatemlwfdist.print(fsulfatemlwfdist);
  watermlwfdist.print(fwatermlwfdist);

  cout << nhso4 << "   " << nso4 << endl;



  fxyz.close();
  fmlwf.close();

}//main

