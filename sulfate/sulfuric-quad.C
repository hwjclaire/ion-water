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
  if (argc<7)
  {
    cout << "USAGE:\n";
    cout << "sulfuric-quad.x [input xyz] [input mlwf] [input quad] #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fquad.open(argv[3]);

  ofstream fdi, fqi, fdw, fqw, fdh, fqh;
  fdi.open("dip-ion");
  fqi.open("quad-ion");
  fdw.open("dip-water");
  fqw.open("quad-water");
  fdh.open("dip-h3o");
  fqh.open("quad-h3o");



  Stat di(0,5,0.01), qi(-10,10,0.01), dw(0,3,0.01);
  Stat qw(-5,4,0.01), dh(0,5,0.01), qh(-6,6,0.01);


  int nhso4 = 0, nso4 = 0;



  int nframe = atoi(argv[4]), natom;
  int nmo = atoi(argv[5]), nskip = atoi(argv[6]);
  Cell c;

  for ( int iframe = 0; iframe < nframe; iframe ++)
  {
    //cout << iframe << endl;

    std::vector<Atom> atomset;
    std::vector<Mlwf> mlwfset;
    std::vector<Water> waterset;


    //read header
    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;
    fxyz >> natom;
    fxyz >> c;

    c.invert();
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
      //cout << atomset[iatom].x()<< endl;

    }//for iatom
  
    //cout << "read\n";

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
    {
      i -> watermove();
      i -> check_atom();
    }

    if ( sulfate.minus1() ) nhso4 ++;
    if ( sulfate.minus2() ) nso4 ++;  

    //cout << "add atom\n";
     
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
        //cout << pwf -> x() << endl;
        //pwf -> readp(fpolar);
        pwf -> readquad(fquad);
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
          min -> add_wf(*pwf);
        }
        else
        {
          sulfate.add_wf(*pwf);
        }
      }//for imo 

      //cout << "add mlwf\n";
      
      sulfate.check_wf();
      sulfate.Compute_cm();
      sulfate.Compute_cc();
      sulfate.Compute_dipole();
      //sulfate.Compute_quad();
      //cout << sulfate.tot_charge() << "  " << sulfate.cm() << "  " << sulfate.dipole() << endl;

      di.add ( length ( sulfate.dipole() ) );
      for ( int k = 0; k < 3 ; k ++ )
        qi.add ( sulfate.quad_eigval()[k] );

  
      for ( vector<Water>::iterator i = waterset.begin(); i < waterset.end(); i ++ )
      {
        i -> mlwfmove();
        i -> check_wf();
        i -> Compute_cm();
        i -> Compute_cc();
        i -> Compute_dipole();
        //i -> Compute_quad();

        //cout << i -> o() -> x() << " " << i -> h(0) -> x() << " " << i -> h(1) -> x() << endl;
        //cout << i->tot_charge() << "  " << i->cm() << "  " << i->dipole() << endl;
        //cout << i -> dipole() << endl;
        if ( i -> isoh3() )
        {
          dh.add ( length( i -> dipole() ) );
          for ( int k = 0; k < 3 ; k ++ )
            qh.add ( i->quad_eigval()[k] );
        }
        else
        {
          dw.add ( length( i -> dipole() ) );
          for ( int k = 0; k < 3 ; k ++ )
            qw.add ( i->quad_eigval()[k] );
        }//if
  
      }//for i

      //cout << "dipole\n";

    }// if iframe % nskip

  }// for iframe

  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  di.print(fdi);
  qi.print(fqi);
  dw.print(fdw);
  qw.print(fqw);
  dh.print(fdh);
  qh.print(fqh);

  cout << nhso4 << "   " << nso4 << endl;



  fxyz.close();
  fmlwf.close();

}//main

