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

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  tm.start();
  if (argc<8)
  {
    cout << "USAGE:\n";
    cout << "dipolemlwf.x [input xyz] [input mlwf] [input polariz]  #frame nMO nskip nskip2\n";
    return 1;
  }

  int nframe = atoi(argv[4]), nmo = atoi(argv[5]), nskip = atoi(argv[6]);
  int nskip2 = atoi(argv[7]), natom;

  //cell
  double volume;
  Cell c;
 
  //input
  ifstream fxyz, fmlwf, fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fpolar.open(argv[3]);

  //output
  ofstream fwfoh, fwflo;
  fwfoh.open("polar_oh");
  fwflo.open("polar_lo");
  fwflo.precision(8);
  fwfoh.precision(8);

  ofstream fpolarstat, fpolarstatavg;
  fpolarstat.open("polarstat");
  fpolarstatavg.open("polarstatavg");

  ofstream fpolarohstat, fpolarlostat;
  fpolarlostat.open("polarohstat");
  fpolarohstat.open("polarlostat");
  
  ofstream polarmol[nmo/4];
  ofstream polaroh[nmo/4];
  ofstream polarlo[nmo/4];

  //statistics
  Stat polarstat(7,23);
  Stat polarstatavg(7,23);
  Stat polarohstat(0,15);
  Stat polarlostat(0,15);

  //data structures
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;


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
        if ( name.compare("O") == 0 ) 
        {
          waterset.push_back( Water(icount ++, c ) );
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
        }
        else
        {
          atomset[iatom].setmass(2.014);
          atomset[iatom].setname(name);
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

      for ( int i = 0; i < watercount; i ++ )
      {
        ostringstream name;
        name << "polarwater" << i+1 << ".dat";  
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmol[i].open(fptr);
        polarmol[i].precision(8);
      
        name.str("");
        name.clear();
        name << "polaroh" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        polaroh[i].open(fptr);
        polaroh[i].precision(8);

        name.str("");
        name.clear();
        name << "polarlo" << i+1 << ".dat";
        filename = name.str();
        fptr = (char*)filename.c_str();
        polarlo[i].open(fptr);
        polarlo[i].precision(8);

      }
    }// if ifram == 0;  
     
    if ( nskip == 0 || iframe % nskip == 0 )
    {
      for ( int i = 0; i < waterset.size(); i ++ )
      {  
        if ( iframe % ( nskip * nskip2 ) != 0 ) continue;
        waterset[i].watermove();
        waterset[i].reset_wf();
      }


      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
        pwf -> readx(fmlwf);
        pwf -> readp(fpolar);
        //cout << imo << " " << pwf -> x() << endl;

        if ( iframe % ( nskip * nskip2 ) != 0 ) continue;
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
        min -> add_wf(*pwf,false);
        //cout << imo << " " << min_dist << endl;
      }
  
      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].check_wf();
      }

#if 1
      if ( iframe % ( nskip * nskip2 ) != 0 ) continue;
      //output polarizabilities
      double sum[6], sumoh[6], sumlo[6];//wf of OH bond and wf of lone pair
      double summol[6];//wf of water molecules.
      double summollo[6];//wf of water molecules lone pair.
      double summoloh[6];//wf of water molecules oh bond.
      //set zero;
      for ( int i = 0; i < 6; i ++)
      {
        sumoh[i] = 0;
        sumlo[i] = 0;
      }

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        const double* p0 = waterset[iwater].wf(0)->polar();
        const double* p1 = waterset[iwater].wf(1)->polar();
        const double* p2 = waterset[iwater].wf(2)->polar();
        const double* p3 = waterset[iwater].wf(3)->polar();

        for ( int i = 0; i < 6; i ++)
        {
          sumoh[i] += p0[i];
          sumoh[i] += p1[i];
          sumlo[i] += p2[i];
          sumlo[i] += p3[i];
          summol[i] = p0[i] + p1[i] + p2[i] + p3[i];
          summoloh[i] = p0[i] + p1[i]; 
          summollo[i] = p2[i] + p3[i]; 


        }

      polarstat.add(summol[0]);
      polarstat.add(summol[3]);
      polarstat.add(summol[5]);
      polarstatavg.add((summol[0]+summol[3]+summol[5])/3);
      

      polarmol [iwater] << summol[0] << "  " << summol[3] << "  " << summol[5] << "  "  
                        << summol[1] << "  " << summol[2] << "  " << summol[4] << "  " 
                        << summol[1] << "  " << summol[2] << "  " << summol[4] << "  " << endl;
      polarlo [iwater]  << summollo[0] << "  " << summollo[3] << "  " << summollo[5] << "  "  
                        << summollo[1] << "  " << summollo[2] << "  " << summollo[4] << "  "
                        << summollo[1] << "  " << summollo[2] << "  " << summollo[4] << "  " << endl;
      polaroh [iwater]  << summoloh[0] << "  " << summoloh[3] << "  " << summoloh[5] << "  "  
                        << summoloh[1] << "  " << summoloh[2] << "  " << summoloh[4] << "  " 
                        << summoloh[1] << "  " << summoloh[2] << "  " << summoloh[4] << "  " << endl;
      }

      //for ( int i = 0; i < 6; i ++) sum[i] = sumoh[i] + sumlo[i];
      fwfoh << sumoh[0] << "  " << sumoh[3] << "  " << sumoh[5] << "  "  
            << sumoh[1] << "  " << sumoh[2] << "  " << sumoh[4] << "  "
            << sumoh[1] << "  " << sumoh[2] << "  " << sumoh[4] << "  " << endl;
      fwflo << sumlo[0] << "  " << sumlo[3] << "  " << sumlo[5] << "  "  
            << sumlo[1] << "  " << sumlo[2] << "  " << sumlo[4] << "  "
            << sumlo[1] << "  " << sumlo[2] << "  " << sumlo[4] << "  " << endl;
#endif


    }// if iframe % nskip

  }// for iframe

  tm.stop();
  cout << "Real Time: " << tm.real() <<  "s  CPU Time:" << tm.cpu() << "s" << endl;

  polarstat.print(fpolarstat);
  polarstatavg.print(fpolarstatavg);

  fxyz.close();
  fmlwf.close();
  fpolar.close();
  fwfoh.close();
  fwflo.close();
  fpolarstat.close();

}//main
