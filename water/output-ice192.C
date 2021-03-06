// output-ice192.x
// For a 192-D2O ice Ih surface (4bilayer) simulation
// output dipole and polarizability of individual water molecules
// in the first bilayer as well.
// Quan (Andy) Wan Thu Nov 14 15:58:13 PST 2013


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
#include "Atom.h"
#include "Mlwf.h"
#include "Timer.h"
#include "Cell.h"
#include "Stat.h"
#include "Tensor.h"

using namespace std;

int main (int argc, char *argv[])
{
  Timer tm;
  Timer tm_addwf, tm_alpha, tm_output, tm_read, tm_misc, tm_rot;
  tm.start();
  if (argc<8)
  {
    cout << "USAGE:\n";
    cout << "output-ice192.x [input xyz] [input mlwf] [input quadrupole] [input polar] #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad,fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  //fquad.open(argv[3]);
  fpolar.open(argv[4]);



  int nframe = atoi(argv[5]), nmo = atoi(argv[6]), nskip = atoi(argv[7]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  ofstream polarmlwf[192];

  ofstream fup,fbo,fup2,fbo2,fmid,fdipole;
  fup.open("surface1");
  fbo.open("surface2");
  fup2.open("subsurface1");
  fbo2.open("subsurface2");
  fmid.open("bulk");
  fdipole.open("dipole");

  Stat quad(-5,5,0.01);
  Stat dip(0,2,0.01);

  ofstream fdip,fq;
  fdip.open("statdip");
  //fq.open("statquad");

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
      //cout << atomset[iatom].x()<< endl;

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


      assert ( waterset.size() * 4 == nmo );
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
      tm_misc.stop();
      //cout << "resetwf\n";

      tm_addwf.start();
      for ( int imo = 0; imo < nmo; imo++)
      {
        Mlwf * pwf = & mlwfset[imo];
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
        min -> add_wf(*pwf,true);
        //cout << imo << " " <<  min_dist << endl;
      }
      tm_addwf.stop();
      //cout << "addwf\n"; 

      tm_misc.start();

      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        waterset[i].mlwfmove();
        waterset[i].check_wf();
        waterset[i].Compute_cc();
        waterset[i].Compute_cm();
        waterset[i].Compute_dipole();
        //waterset[i].Compute_quad();
        waterset[i].Compute_polariz();
        //for ( int k = 0; k < 3 ; k ++ )
        //  quad.add ( waterset[i].quad_eigval()[k] );
        dip.add(length(waterset[i].dipole()));

        //output each molecules
        polarmlwf[i] << waterset[i].dipole();
        waterset[i].polariz().print_inline(polarmlwf[i]);
        polarmlwf[i] << endl;
      }
      tm_misc.stop();
      //cout << "Check Done" << endl;
  
      tm_output.start();

#if 1

      D3vector d1(0,0,0),d2(0,0,0),d3(0,0,0),d4(0,0,0),d5(0,0,0),d6(0,0,0);
      Tensor p1,p2,p3,p4,p5,p6;
      int count1=0,count2=0,count3=0,count4=0,count5=0,count6=0;
     

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        double z = waterset[iwater].o()->x().z;
        //cout << z << endl;
        if ( z < 6 ) 
        {
          d1 += waterset[iwater].dipole();
          p1 += waterset[iwater].polariz();
          count1++;
        }
        else if ( z < 14 )
        {
          d2 += waterset[iwater].dipole();
          p2 += waterset[iwater].polariz();
          count2++;
        }
        else if ( z < 20 )
        {
          d3 += waterset[iwater].dipole();
          p3 += waterset[iwater].polariz();
          count3++;
        }
        else 
        {
          d4 += waterset[iwater].dipole();
          p4 += waterset[iwater].polariz();
          count4++;
        }

      }//for iwater

      assert( count1 == 48 );
      assert( count2 == 48 );
      assert( count3 == 48 );
      assert( count4 == 48 );
      
      d1.z=-d1.z;
      d2.z=-d2.z;

      fbo << d1;
      p1.print_inline(fbo);
      fbo << endl;

      fbo2 << d2;
      p2.print_inline(fbo2);
      fbo2 << endl;

      fup << d4;
      p4.print_inline(fup);
      fup << endl;

      fup2 << d3;
      p3.print_inline(fup2);
      fup2 << endl;
  
      D3vector d0(0,0,0), dmid(0,0,0);
      Tensor p0,pmid;
      d0 = d1+d2+d3+d4;
      p0 = p1+p2+p3+p4;

      dmid=d2-d3;
      pmid=p2-p3;

      fdipole << d0 << endl;
      
      fmid << dmid;
      pmid.print_inline(fmid);
      fmid << endl;


#endif
      tm_output.stop();

    }// if iframe % nskip

  }// for iframe

  //quad.print(fq);
  dip.print(fdip);

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
  //fquad.close();

}//main
