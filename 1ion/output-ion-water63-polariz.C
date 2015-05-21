// output-ion-water63-polariz.x
// For a 63-D2O 1-ion  simulation
// and compute molecular polarizabilities using ewald sum
// Quan (Andy) Wan Mon May 18 18:17:05 CDT 2015


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
#include "Mol.h"
#include "Water.h"
#include "Halide.h"
#include "Alkali.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "Stat.h"
#include "Tensor.h"
#include "Polariz.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<8)
  {
    cout << "USAGE:\n";
    cout << "output-ion-water63-polarize.x [input xyz] [input mlwf] [input quadrupole] [input polar] #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad,fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fquad.open(argv[3]);
  fpolar.open(argv[4]);

  Polariz * polariz = 0;

  int nframe = atoi(argv[5]), nmo = atoi(argv[6]), nskip = atoi(argv[7]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;
  std::vector<Halide> halideset;
  std::vector<Alkali> alkaliset;

  
  ofstream polarmlwf[64];

  Stat alphastat1(5,25,0.01);
  Stat alphastat2(5,25,0.01);
  Stat alphastat3(5,25,0.01);

  ofstream falphastat1, falphastat2, falphastat3;
  falphastat1.open("statalpha1");
  falphastat2.open("statalpha2");
  falphastat3.open("statalpha3");

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
    int wcount = 0, hcount = 0, acount =0;
    tm_read.start();
    for ( int iatom = 0; iatom < natom; iatom ++ )
    {
      string name;

      fxyz >> name;

      if ( iframe == 0 )
      {
        if ( name.compare("Cl") == 0 )
        {
          halideset.push_back( Halide(hcount ++, c ) );
          atomset[iatom].setmass(35.45);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(7.0);
        }
        else if ( name.compare("Na") == 0 )
        {
          alkaliset.push_back( Alkali(acount ++, c ) );
          atomset[iatom].setmass(22.99);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6.0);
        }
        else if ( name.compare("O") == 0 ) 
        {
          waterset.push_back( Water(wcount ++, c ) );
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
      //cout << atomset[iatom].x()<< endl;

    }//for iatom
    tm_read.stop();

    //assign H to O
    if ( iframe == 0 )  //assign molecules 
    {

      cout << "natom " << atomset.size() << endl;
      cout << "nwater " << waterset.size() << endl;


      int watercount = 0;
      int halidecount = 0;
      int alkalicount = 0;
      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        if ( atomset[iatom].name().compare("Cl") == 0 )  // is Cl
        {
          halideset[halidecount].add_halide(atomset[iatom]);
          halidecount ++;
          continue;
        }

        if ( atomset[iatom].name().compare("Na") == 0 )  // is Na
        {
          alkaliset[alkalicount].add_alkali(atomset[iatom]);
          alkalicount ++;
          continue;
        }

        if ( atomset[iatom].name().compare("O")  ) continue;//is not O
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

      //assert ( halidecount + alkalicount + watercount == 64 );
      cout << "Number of water  " << watercount << endl;
      cout << "Number of halide " << halidecount << endl;
      cout << "Number of alkali " << alkalicount << endl;


      assert ( waterset.size() * 4 + halideset.size() * 4 
               + alkaliset.size() * 3 == nmo );
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        waterset[i].check_atom();
      }
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        ostringstream name;
        name << "water" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmlwf[i].open(fptr);
        polarmlwf[i].setf(ios::fixed, ios::floatfield);
        polarmlwf[i].setf(ios::right, ios::adjustfield);
        polarmlwf[i].precision(5);
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
        //pwf -> readquad(fquad);
        pwf -> readp(fpolar);
        pwf -> setnumber(imo);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();



      tm_misc.start();
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        //update oh distances.
        waterset[i].watermove();
        waterset[i].reset_wf();
      }
      for ( int i = 0; i < halideset.size(); i ++ )
      {
        halideset[i].reset_wf();
      }
      for ( int i = 0; i < alkaliset.size(); i ++ )
      {
        alkaliset[i].reset_wf();
      }




      tm_misc.stop();
      //cout << "resetwf\n";


      tm_addwf.start();
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];

        // distance to a water
        double min_dist = 100;
        Water * min;
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

        // distance to a halide
        double  hmin_dist = 100;
        Halide * hmin;
        for ( int ihalide = 0; ihalide < halideset.size(); ihalide ++)
        {
          Halide * phalide = &halideset[ihalide];
          //if (waterset[iwater].wf_full()) continue;
          double dist = phalide -> distance(pwf);
          if ( dist < hmin_dist )
          {
            hmin = phalide;
            hmin_dist = dist;
          }
        }

        // distance to a alkali
        double  amin_dist = 100;
        Alkali * amin;
        for ( int ialkali = 0; ialkali < alkaliset.size(); ialkali ++)
        {
          Alkali * palkali = &alkaliset[ialkali];
          //if (waterset[iwater].wf_full()) continue;
          double dist = palkali -> distance(pwf);
          if ( dist < hmin_dist )
          {
            amin = palkali;
            amin_dist = dist;
          }
        }


        if ( min_dist < hmin_dist && min_dist < amin_dist )
        {
          min -> add_wf(*pwf,true);
        }
        else if ( hmin_dist < min_dist && hmin_dist < amin_dist )
        {
          hmin -> add_wf(*pwf,true);
        }
        else
        {
          amin -> add_wf(*pwf,true);
        }
      }
      tm_addwf.stop();
      //cout << "addwf\n"; 

      tm_misc.start();


      // compute polarizabilities for water
      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].mlwfmove();
        waterset[i].check_wf();
        waterset[i].Compute_cc();
        waterset[i].Compute_cm();
        waterset[i].Compute_dipole();
        //waterset[i].Compute_quad();
        waterset[i].Compute_polariz();
  
        //output each molecules
        //polarmlwf[i] << waterset[i].dipole();
        //waterset[i].polariz().print_inline(polarmlwf[i]);
        //polarmlwf[i] << endl;
      }

      // compute polarizabilities of halide
      for ( int i = 0; i < halideset.size(); i ++ )
      {
        halideset[i].check_wf();
        halideset[i].Compute_cc();
        halideset[i].Compute_cm();
        halideset[i].Compute_polariz();
      }

      // compute polarizabilities of alkali
      for ( int i = 0; i < alkaliset.size(); i ++ )
      {
        alkaliset[i].check_wf();
        alkaliset[i].Compute_cc();
        alkaliset[i].Compute_cm();
        alkaliset[i].Compute_polariz();
      }



      tm_misc.stop();
      //cout << "Check Done" << endl;
 
      if ( iframe == 0 )
      {

        // copy waterset and halideset and alkaliset to molset for ewald
        molset.resize(waterset.size()+halideset.size()+alkaliset.size());
        for ( int i = 0; i < waterset.size(); i ++ )
        {
          molset[i] = (Mol*)&waterset[i];
          //cout << molset[i] -> cc() << endl;
        }

        for ( int i = waterset.size(); i < waterset.size()+halideset.size(); i ++ )
        {
          molset[i] = (Mol*)&halideset[i-waterset.size()];
          //cout << molset[i] -> cc() << endl;
        }


        for ( int i = waterset.size()+halideset.size(); i < molset.size(); i ++ )
        {
          molset[i] = (Mol*)&alkaliset[i-waterset.size()-halideset.size()];
          //cout << molset[i] -> cc() << endl;
        }


        polariz = new Polariz(c, molset, mlwfset);
      }
      //cout << "Init Ewald" << endl;

      //cout << molset.size() << endl;

      tm_alpha.start();
      polariz -> Calculate_mol_Ewald();
      tm_alpha.stop();
      //cout << "Ewald done" << endl;

      tm_output.start();
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        waterset[i].getaxis();
        waterset[i].rotate_alpha_tensor();
        waterset[i].rotate_polar_tensor();
        alphastat1.add(waterset[i].alpha_tensor_axis(0));
        alphastat2.add(waterset[i].alpha_tensor_axis(4));
        alphastat3.add(waterset[i].alpha_tensor_axis(8));
      }


      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        polarmlwf[iwater] 
          << setw(12) << waterset[iwater].alpha_tensor_axis()[0]
          << setw(12) << waterset[iwater].alpha_tensor_axis()[4]
          << setw(12) << waterset[iwater].alpha_tensor_axis()[8]
          << setw(12) << waterset[iwater].alpha_tensor_axis()[1]
          << setw(12) << waterset[iwater].alpha_tensor_axis()[2]
          << setw(12) << waterset[iwater].alpha_tensor_axis()[5]
          << endl;
      }



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
  fquad.close();

  cout << alphastat1.n() << "   " << alphastat1.avg() << endl;
  cout << alphastat2.n() << "   " << alphastat2.avg() << endl;
  cout << alphastat3.n() << "   " << alphastat3.avg() << endl;

  alphastat1.print(falphastat1);
  alphastat2.print(falphastat2);
  alphastat3.print(falphastat3);



}//main
