#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "D3vector.h"
#include "Mol.h"
#include "NG.h"
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "Stat.h"
#include "Polariz.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  tm.start();
  if (argc<7)
  {
    cout << "USAGE:\n";
    cout << "testPolariz.x [input xyz] [input mlwf] [input polariz]  #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fpolar.open(argv[3]);


  ofstream alphamolaxis[2];
  ofstream alphamol[2];
  ofstream polarmol[2];
  
  ofstream polarmlwfaxis[128];

  ofstream dipole;
  dipole.open("totaldipole.dat");

  int nframe = atoi(argv[4]), nmo = atoi(argv[5]), nskip = atoi(argv[6]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<NG> NGset;
  std::vector<Mol*> molset;

  Polariz * polariz = 0;

  for ( int iframe = 0; iframe < nframe; iframe ++)
  {

    //read header
    if ( iframe % 1000 == 0 ) cout << "Frame" << iframe << endl;
    fxyz >> natom;
    fxyz >> c;
  
    if ( iframe == 0 )
    {
      c.invert();
      atomset.resize(natom);
      mlwfset.resize(nmo);
    }

    //read atoms
    int icount = 0;
    for ( int iatom = 0; iatom < natom; iatom ++ )
    {
      string name;
      fxyz >> name;
      if ( iframe == 0 )
      {
        if ( name.compare("Ar") == 0 || name.compare("Ne") == 0) 
        {
          NGset.push_back( NG(icount ++, c ) );

          if ( name.compare("Ar") == 0)
            atomset[iatom].setmass(39.948);
          else
            atomset[iatom].setmass(20.1797);

          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6.0);
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].pos()<< endl;

    }//for iatom

    //assign H to O
    if ( iframe == 0 )  //assign molecules 
    {

      cout << "natom " << atomset.size() << endl;
      cout << "nNG " << NGset.size() << endl;

      int ngcount = 0;
      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        //cout << NGcount << & NGset[NGcount] << " " << & atomset[iatom] << endl;
        NGset[ngcount].add_oxygen ( atomset[iatom] );
        ngcount ++;
      }

      for ( int i = 0; i < NGset.size(); i ++ )
      {
        NGset[i].check_atom();
      }


      for ( int i = 0; i < ngcount; i ++ )
      {
        ostringstream name;
        name << "alphaNG" << i+1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        alphamol[i].open(fptr);
        alphamol[i].setf(ios::fixed, ios::floatfield);
        alphamol[i].setf(ios::right, ios::adjustfield);
        alphamol[i].precision(8);

        name.str("");
        name.clear();
        name << "alphaNGaxis" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        alphamolaxis[i].open(fptr);
        alphamolaxis[i].setf(ios::fixed, ios::floatfield);
        alphamolaxis[i].setf(ios::right, ios::adjustfield);
        alphamolaxis[i].precision(8);




      }
    }// if ifram == 0;  
    //cout << "addh\n";

    if ( nskip == 0 || iframe % nskip == 0 )
    {
      
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x() << endl;
      }
      if ( iframe % ( nskip * 4 ) != 0 ) continue;

      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        double min_dist = 100;
        Mol * min;
        int nNG = NGset.size();
        for ( int iNG = 0; iNG < nNG; iNG ++)
        {
          NG * pNG = &NGset[iNG];
          //if (NGset[iNG].wf_full()) continue;
          double dist = pNG -> distance(pwf);
          if ( dist < min_dist )
          {
            min = pNG;
            min_dist = dist;
          }
        }
        min -> add_wf(*pwf);
        //cout << imo << " " << min_dist << endl;
      }

      for ( int i = 0; i < NGset.size(); i ++ ) 
      {
        NGset[i].check_wf();
        NGset[i].Compute_cm();
        NGset[i].Compute_cc();
        NGset[i].Compute_dipole();
      }
      //cout << "Check Done" << endl;

      if ( iframe == 0 )
      {
        molset.resize(NGset.size());
        for ( int i = 0; i < molset.size(); i ++ )
          molset[i] = (Mol*)&NGset[i];

        polariz = new Polariz(c, molset, mlwfset);
      }
    
      //cout << "polar" << endl;
      //polariz -> Calculate_water();
      polariz -> Calculate_mol_Ewald();
      //polariz -> Calculate_mlwf_Ewald();
      //polariz -> Polar_mlwf_rotate();
      //polariz -> Polar_water_rotate();

#if 0

      for ( int iNG = 0; iNG < NGset.size(); iNG ++ )
      for ( int imlwf = 0; imlwf < 4; imlwf ++ )
      {
        Mlwf* pwf = NGset[iNG].wf(imlwf);
        polarmlwfaxis [iNG * 4 + imlwf] 
                          << setw(12) << pwf->polar_tensor_axis()[0]
                          << setw(12) << pwf->polar_tensor_axis()[4]
                          << setw(12) << pwf->polar_tensor_axis()[8]
                          << setw(12) << pwf->polar_tensor_axis()[1]
                          << setw(12) << pwf->polar_tensor_axis()[2]
                          << setw(12) << pwf->polar_tensor_axis()[5]
                          << setw(12) << pwf->polar_tensor_axis()[3]
                          << setw(12) << pwf->polar_tensor_axis()[6]
                          << setw(12) << pwf->polar_tensor_axis()[7]
                          << endl;
        //cout << iNG * 4 + imlwf << endl;

        for ( int i = 0; i < 9; i ++ ) 
        {
          //cout << pwf->polar_tensor_axis()[i] << endl;
          if ( imlwf < 2 )
            polarmlwfohaxisstat[i]->add(pwf->polar_tensor_axis()[i]);
          else
            polarmlwfloaxisstat[i]->add(pwf->polar_tensor_axis()[i]);
        }

      }

#endif

#if 0
/*
      for ( int i = 0; i < NGset.size(); i ++ )
      {
        alphastat[0]->add ( NGset[i].polar_tensor_axis(0) );
        alphastat[1]->add ( NGset[i].polar_tensor_axis(4) );
        alphastat[2]->add ( NGset[i].polar_tensor_axis(8) );
      }
*/
      for ( int iNG = 0; iNG < NGset.size(); iNG ++ )
      {
        polarmol [iNG]
                          << setw(12) << NGset[iNG].polar_tensor()[0]
                          << setw(12) << NGset[iNG].polar_tensor()[4]
                          << setw(12) << NGset[iNG].polar_tensor()[8]
                          << setw(12) << NGset[iNG].polar_tensor()[1]
                          << setw(12) << NGset[iNG].polar_tensor()[2]
                          << setw(12) << NGset[iNG].polar_tensor()[5]
                          << setw(12) << NGset[iNG].polar_tensor()[3]
                          << setw(12) << NGset[iNG].polar_tensor()[6]
                          << setw(12) << NGset[iNG].polar_tensor()[7]
                          << endl;
      }

#endif


#if 1
  D3vector tot_dipole(0,0,0);
  for ( int i = 0; i < NGset.size(); i ++ )
  {
    tot_dipole += NGset[i].dipole();
/*
    cout << i + 1 << endl;
    cout << NGset[i].cm() << endl;
    cout << NGset[i].cc() << endl;
    cout << NGset[i].dipole() << endl;
*/
  }

  dipole << tot_dipole << endl;
#endif


#if 1

      for ( int iNG = 0; iNG < NGset.size(); iNG ++ )
      {
        alphamol [iNG] << setw(12) << NGset[iNG].alpha_tensor()[0]
                          << setw(12) << NGset[iNG].alpha_tensor()[4]
                          << setw(12) << NGset[iNG].alpha_tensor()[8]
                          << setw(12) << NGset[iNG].alpha_tensor()[1]
                          << setw(12) << NGset[iNG].alpha_tensor()[2]
                          << setw(12) << NGset[iNG].alpha_tensor()[5]
                          << setw(12) << NGset[iNG].alpha_tensor()[3]
                          << setw(12) << NGset[iNG].alpha_tensor()[6]
                          << setw(12) << NGset[iNG].alpha_tensor()[7]
                          << endl;
      }

#endif

    }// if iframe % nskip

  }// for iframe


  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

#if 1
  for ( int i = 0; i < 3; i++ )
  {
    //cout << alphastat[i]->n() << "   " << alphastat[i]->avg() << endl;

    //alphastat[i]->print(falphastat[i]);  
  }
#endif





#if 0
  for ( int i = 0; i < 9; i++ )
  {
    cout << polarmlwfohaxisstat[i]->n() << "   " << polarmlwfohaxisstat[i]->avg() <<  "   ";
    cout << polarmlwfloaxisstat[i]->n() << "   " << polarmlwfloaxisstat[i]->avg() << endl;

    polarmlwfohaxisstat[i]->print(fpolarmlwfohaxisstat[i]);  
    polarmlwfloaxisstat[i]->print(fpolarmlwfloaxisstat[i]);  
  }
#endif

  fxyz.close();
  fmlwf.close();
  fpolar.close();

}//main
