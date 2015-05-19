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
#include "Water.h"
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
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<6)
  {
    cout << "USAGE:\n";
    cout << "PolarizDM.x [input xyz] [input mlwf] [NOT used]  #frame nMO\n";
    return 1;
  }

  int nframe = atoi(argv[4]), nmo = atoi(argv[5]);
  double volume;
  Cell c;
 
  ifstream fxyz, fmlwf, fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fpolar.open(argv[3]);

  

  Stat * alphastat[3];
  for ( int i = 0; i < 3; i ++ ) alphastat[i] = new Stat(5,25);
  ofstream falphastat[3]; 
  falphastat[0].open("alphastat1.dat");
  falphastat[1].open("alphastat2.dat");
  falphastat[2].open("alphastat3.dat");


  Stat alphastatavg(7,22);

  Stat* polarmlwfohaxisstat[9], *polarmlwfloaxisstat[9];

  ofstream fpolarmlwfohaxisstat[9], fpolarmlwfloaxisstat[9];
  for ( int i = 0; i < 9; i ++ )
  {
    polarmlwfohaxisstat[i] = new Stat(-15,15);
    polarmlwfloaxisstat[i] = new Stat(-15,15);


    ostringstream name;
    name << "polarmlwfohaxisstat" << i+1 << ".dat";
    string filename = name.str();
    char * fptr = (char*)filename.c_str();
    //fpolarmlwfohaxisstat[i].open(fptr);

    name.str("");
    name.clear();
    name << "polarmlwfloaxisstat" << i+1 << ".dat";
    filename = name.str();
    fptr = (char*)filename.c_str();
    //fpolarmlwfloaxisstat[i].open(fptr);
  }

  ofstream polarmlwfaxis[mlwf];

  std::vector<Mlwf> mlwfset;

  Polariz * polariz = 0;

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

    //cout << "Read atom\n";

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
      //cout << "addh\n";

      for ( int i = 0; i < nmo; i ++ )
      {
        ostringstream name;
        name << "polarmlwf" << i+1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        //polarmlwfaxis[i].open(fptr);
        //polarmlwfaxis[i].setf(ios::fixed, ios::floatfield);
        //polarmlwfaxis[i].setf(ios::right, ios::adjustfield);
        //polarmlwfaxis[i].precision(8);
      }
      //cout << "addh\n";


      for ( int i = 0; i < watercount; i ++ )
      {
        ostringstream name;
        name << "polarwateraxis" << i+1 << ".dat";  
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmolaxis[i].open(fptr);
        polarmolaxis[i].setf(ios::fixed, ios::floatfield);
        polarmolaxis[i].setf(ios::right, ios::adjustfield);

        name.str("");
        name.clear();
        name << "polarwater" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        polarmol[i].open(fptr);
        polarmol[i].setf(ios::fixed, ios::floatfield);
        polarmol[i].setf(ios::right, ios::adjustfield);
        polarmol[i].precision(4);

        name.str("");
        name.clear();
        name << "dipolewater" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        dipolemol[i].open(fptr);
        dipolemol[i].setf(ios::fixed, ios::floatfield);
        dipolemol[i].setf(ios::right, ios::adjustfield);
        dipolemol[i].precision(4);

        name.str("");
        name.clear();
        name << "alphawater" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        alphamol[i].open(fptr);
        alphamol[i].setf(ios::fixed, ios::floatfield);
        alphamol[i].setf(ios::right, ios::adjustfield);
        alphamol[i].precision(4);

        name.str("");
        name.clear();
        name << "alphawateraxis" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        alphamolaxis[i].open(fptr);
        alphamolaxis[i].setf(ios::fixed, ios::floatfield);
        alphamolaxis[i].setf(ios::right, ios::adjustfield);
        alphamolaxis[i].precision(4);

      


        name.str("");
        name.clear();
        name << "polaroh" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        polaroh[i].open(fptr);
        polaroh[i].setf(ios::fixed, ios::floatfield);
        polaroh[i].setf(ios::right, ios::adjustfield);
        polaroh[i].precision(4);

  

        name.str("");
        name.clear();
        name << "polarlo" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        polarlo[i].open(fptr);
        polarlo[i].setf(ios::fixed, ios::floatfield);
        polarlo[i].setf(ios::right, ios::adjustfield);
        polarlo[i].precision(4);


      }
    }// if ifram == 0;  

    if ( nskip == 0 || iframe % nskip == 0 )
    {

      tm_read.start();     
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();

      if ( iframe % ( nskip * nskip2 ) != 0 ) continue;
  

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
        waterset[i].Compute_cm();
        //waterset[i].Compute_cc();
        waterset[i].getaxis();
        waterset[i].Compute_dipole();
      }
      tm_misc.stop();
      //cout << "Check Done" << endl;

      if ( iframe == 0 )
      {
        molset.resize(waterset.size());
        for ( int i = 0; i < molset.size(); i ++ )
          molset[i] = (Mol*)&waterset[i];

        polariz = new Polariz(c, molset, mlwfset);
      }
    
      //cout << "polar" << endl;
      //polariz -> Calculate_water();
      tm_alpha.start();
      polariz -> Calculate_mol_Ewald();
      tm_alpha.stop();
      //polariz -> Calculate_mlwf_Ewald();
      //polariz -> Polar_mlwf_rotate();
      //polariz -> Polar_water_rotate();

      //cout << " Ewald Done" << endl;

      tm_rot.start();
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        waterset[i].rotate_alpha_tensor();
        waterset[i].rotate_polar_tensor();
      }
      tm_rot.stop();


      tm_output.start();
#if 0

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      for ( int imlwf = 0; imlwf < 4; imlwf ++ )
      {
        Mlwf* pwf = waterset[iwater].wf(imlwf);
        polarmlwfaxis [iwater * 4 + imlwf] 
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
        //cout << iwater * 4 + imlwf << endl;

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


#if 1
      std::vector<double> sumlo;
      std::vector<double> sumoh;
      sumlo.resize(6,0);
      sumoh.resize(6,0);

      for ( int i = 0; i < waterset.size(); i ++ )
      {
        double summoloh[6];
        double summollo[6];
        for ( int j = 0; j < 6; j ++ )
        {
          summoloh[j] = waterset[i].wf(0)->polar()[j] + waterset[i].wf(1)->polar()[j];
          summollo[j] = waterset[i].wf(2)->polar()[j] + waterset[i].wf(3)->polar()[j];

          sumoh[j] += summoloh[j];
          sumlo[j] += summollo[j];

        }

        polarlo [i]  << summollo[0] << "  " << summollo[3] << "  " << summollo[5] << "  "  
                     << summollo[1] << "  " << summollo[2] << "  " << summollo[4] << "  "
                     << summollo[1] << "  " << summollo[2] << "  " << summollo[4] << "  " << endl;
        polaroh [i]  << summoloh[0] << "  " << summoloh[3] << "  " << summoloh[5] << "  "  
                     << summoloh[1] << "  " << summoloh[2] << "  " << summoloh[4] << "  " 
                     << summoloh[1] << "  " << summoloh[2] << "  " << summoloh[4] << "  " << endl;


      }

      fwfoh << sumoh[0] << "  " << sumoh[3] << "  " << sumoh[5] << "  "
            << sumoh[1] << "  " << sumoh[2] << "  " << sumoh[4] << "  "
            << sumoh[1] << "  " << sumoh[2] << "  " << sumoh[4] << "  " << endl;
      fwflo << sumlo[0] << "  " << sumlo[3] << "  " << sumlo[5] << "  "
            << sumlo[1] << "  " << sumlo[2] << "  " << sumlo[4] << "  "
            << sumlo[1] << "  " << sumlo[2] << "  " << sumlo[4] << "  " << endl;

#endif
#if 1
/*
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        alphastat[0]->add ( waterset[i].polar_tensor_axis(0) );
        alphastat[1]->add ( waterset[i].polar_tensor_axis(4) );
        alphastat[2]->add ( waterset[i].polar_tensor_axis(8) );
      }
*/

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        polarmol [iwater]
                          << setw(12) << waterset[iwater].polar_tensor()[0]
                          << setw(12) << waterset[iwater].polar_tensor()[4]
                          << setw(12) << waterset[iwater].polar_tensor()[8]
                          << setw(12) << waterset[iwater].polar_tensor()[1]
                          << setw(12) << waterset[iwater].polar_tensor()[2]
                          << setw(12) << waterset[iwater].polar_tensor()[5]
                          << setw(12) << waterset[iwater].polar_tensor()[3]
                          << setw(12) << waterset[iwater].polar_tensor()[6]
                          << setw(12) << waterset[iwater].polar_tensor()[7]
                          << endl;
      }



      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        polarmolaxis [iwater]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[0]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[4]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[8]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[1]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[2]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[5]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[3]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[6]
                          << setw(12) << waterset[iwater].polar_tensor_axis()[7]
                          << endl;
      }

#endif


#if 1
  D3vector tot_dipole(0,0,0);
  for ( int i = 0; i < waterset.size(); i ++ )
  {

    dipolemol[i] << waterset[i].dipole() << endl;
    tot_dipole += waterset[i].dipole();
/*
    cout << i + 1 << endl;
    cout << waterset[i].cm() << endl;
    cout << waterset[i].cc() << endl;
    cout << waterset[i].dipole() << endl;
*/
  }

  dipole << tot_dipole << endl;
#endif


#if 1

      for ( int i = 0; i < waterset.size(); i ++ )
      {
        alphastat[0]->add ( waterset[i].alpha_tensor_axis(0) );
        alphastat[1]->add ( waterset[i].alpha_tensor_axis(4) );
        alphastat[2]->add ( waterset[i].alpha_tensor_axis(8) );
      }
  

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        alphamol [iwater] << setw(12) << waterset[iwater].alpha_tensor()[0]
                          << setw(12) << waterset[iwater].alpha_tensor()[4]
                          << setw(12) << waterset[iwater].alpha_tensor()[8]
                          << setw(12) << waterset[iwater].alpha_tensor()[1]
                          << setw(12) << waterset[iwater].alpha_tensor()[2]
                          << setw(12) << waterset[iwater].alpha_tensor()[5]
                          << setw(12) << waterset[iwater].alpha_tensor()[3]
                          << setw(12) << waterset[iwater].alpha_tensor()[6]
                          << setw(12) << waterset[iwater].alpha_tensor()[7]
                          << endl;
      }

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        alphamolaxis [iwater] 
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[0]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[4]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[8]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[1]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[2]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[5]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[3]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[6]
                          << setw(12) << waterset[iwater].alpha_tensor_axis()[7]
                          << endl;
      }

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

#if 1
  for ( int i = 0; i < 3; i++ )
  {
    cout << alphastat[i]->n() << "   " << alphastat[i]->avg() << endl;

    alphastat[i]->print(falphastat[i]);  
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
  fwfoh.close();
  fwflo.close();

}//main
