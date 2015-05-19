// output-ice144-quad-list.x
// For a 144-D2O ice Ih surface simulation ( 3 bilayers )
// added output for list of water molecules in each surfaces
// Quan (Andy) Wan Thu Jan  1 17:46:21 CST 2015


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
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
    cout << "output-ice.x [input xyz] [input mlwf] [input quadrupole] [input polar] #frame nMO nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad,fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fquad.open(argv[3]);
  fpolar.open(argv[4]);



  int nframe = atoi(argv[5]), nmo = atoi(argv[6]), nskip = atoi(argv[7]), natom;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  ofstream polarmlwf[144];

  ofstream fup,fbo,fup2,fbo2,fup3,fbo3,fmid,fdipole,foxyz;
  fup.open("surface1");
  fbo.open("surface2");
  fup2.open("subsurface1");
  fbo2.open("subsurface2");
  fup3.open("ssurface1");
  fbo3.open("ssurface2");
  fmid.open("bulk");  
  fdipole.open("dipole");
  //foxyz.open("out.xyz");

  //Stat quad(-5,5,0.01);
  //Stat pol(5,18,0.01);
  //Stat dip(0,2,0.01);
  //Stat owf(0,2,0.01);

  ofstream fdip,fq,fowf,fpol;
  //fdip.open("statdip");
  //fpol.open("statpol");
  //fq.open("statquad");
  //fowf.open("statowf");

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
        pwf -> readquad(fquad);
        pwf -> readp(fpolar);
        pwf -> setnumber(imo);
        //cout << imo << " " << pwf -> x() << endl;
      }
      tm_read.stop();

#if 0 
      // output xyz including atoms and MLWFs
      foxyz << natom + nmo << endl;
      for ( int i = 0; i < 9; i ++ )
        foxyz << c.a()[i] << " ";
      foxyz << endl;

      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        D3vector vec = atomset[iatom].x();
        c.images(vec);
        foxyz << atomset[iatom].name() << " "<< vec*.529 << endl;
      }


      for ( int imo = 0; imo < nmo; imo++)
        foxyz << "M " << mlwfset[imo].x()*.529 << endl;

#endif

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
        //owf.add(min_dist);
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
        waterset[i].Compute_quad();
        waterset[i].Compute_polariz();
        for ( int k = 0; k < 3 ; k ++ )
        {
          //quad.add ( waterset[i].quad_eigval()[k] );
          //pol.add( waterset[i].polariz()[4*k] );
        }
        //dip.add(length(waterset[i].dipole()));
      }

      tm_misc.stop();
      //cout << "Check Done" << endl;
  
      tm_output.start();

#if 1
      D3vector d1(0,0,0),d2(0,0,0),d3(0,0,0),d4(0,0,0),d5(0,0,0),d6(0,0,0);
      Tensor p1,p2,p3,p4,p5,p6;
      int count1=0,count2=0,count3=0,count4=0,count5=0,count6=0;
      
      // set of water numbers
      std::vector<int> wl1, wl2, wl3, wl4, wl5, wl6;


      std::map <double,int> zmap;
      std::map <double,int>::iterator imap;
      // sort the z coordinate of all water in waterset
      // zmap.first store z coordinate, zmap.second store index of 
      // the water in waterset
      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {
        double z = waterset[iwater].o()->x().z;
        zmap[z] = iwater;
      }

      // group waters into bilayers based on their z coordinate
      int icount = 0;
      for ( imap = zmap.begin(); imap != zmap.end(); imap ++, icount ++)
      {
        int iwater = imap -> second;
        //cout << imap -> first << " " << iwater << endl;
        if ( icount < 12 ) 
        {
          d1 += waterset[iwater].dipole();
          p1 += waterset[iwater].polariz();
          if ( iframe == 0 ) wl1.push_back(waterset[iwater].number());
          count1++;
          waterset[iwater].dipole().z = waterset[iwater].dipole().z;
        }
        else if ( icount < 48 )
        {
          d2 += waterset[iwater].dipole();
          p2 += waterset[iwater].polariz();
          if ( iframe == 0 ) wl2.push_back(waterset[iwater].number());
          count2++;
          waterset[iwater].dipole().z = waterset[iwater].dipole().z;
        }
        else if ( icount < 72 )
        {
          d3 += waterset[iwater].dipole();
          p3 += waterset[iwater].polariz();
          if ( iframe == 0 ) wl3.push_back(waterset[iwater].number());
          count3++;
          waterset[iwater].dipole().z = waterset[iwater].dipole().z;
        }
        else if ( icount < 96 ) 
        {
          d4 += waterset[iwater].dipole();
          p4 += waterset[iwater].polariz();
          if ( iframe == 0 ) wl4.push_back(waterset[iwater].number());
          count4++;
        }
        else if ( icount < 132 )
        {
          d5 += waterset[iwater].dipole();
          p5 += waterset[iwater].polariz();
          if ( iframe == 0 ) wl5.push_back(waterset[iwater].number());
          count5++;
        }
        else
        {
          d6 += waterset[iwater].dipole();
          p6 += waterset[iwater].polariz();
          if ( iframe == 0 ) wl6.push_back(waterset[iwater].number());
          count6++;
        }

      }
      //cout << count1 << endl;
      //cout << count2 << endl;
      //cout << count3 << endl;
      //cout << count4 << endl;
  
      d1.z = -d1.z;
      d2.z = -d2.z;
      d3.z = -d3.z;


      assert( count1 == 12 );
      assert( count2 == 36 );
      assert( count3 == 24 );
      assert( count4 == 24 );
      assert( count5 == 36 );
      assert( count6 == 12 );


      // couput name list
      if ( iframe == 0 )
      {
        fbo << count1 << endl;
        fbo2 << count2 << endl;
        fbo3 << count3 << endl;
        fup << count6 << endl;
        fup2 << count5 << endl;
        fup3 << count4 << endl;

        std::vector<int>::iterator iwl;
        for ( iwl = wl1.begin(); iwl < wl1.end(); iwl ++ )
          fbo << *iwl << endl;
        for ( iwl = wl2.begin(); iwl < wl2.end(); iwl ++ )
          fbo2 << *iwl << endl;
        for ( iwl = wl3.begin(); iwl < wl3.end(); iwl ++ )
          fbo3 << *iwl << endl;
        for ( iwl = wl6.begin(); iwl < wl6.end(); iwl ++ )
          fup << *iwl << endl;
        for ( iwl = wl5.begin(); iwl < wl5.end(); iwl ++ )
          fup2 << *iwl << endl;
        for ( iwl = wl4.begin(); iwl < wl4.end(); iwl ++ )
          fup3 << *iwl << endl;

      }


      // output water later
      // need to flip dipole first
      for ( int i = 0; i < waterset.size(); i ++ ) 
      {
        //output each molecules
        polarmlwf[i] << waterset[i].cm();
        polarmlwf[i] << waterset[i].dipole();
        waterset[i].quad().print_inline(polarmlwf[i]);
        polarmlwf[i] << endl;
      }

#endif
      tm_output.stop();

    }// if iframe % nskip

  }// for iframe

  //quad.print(fq);
  //dip.print(fdip);
  //pol.print(fpol);
  //owf.print(fowf);

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

}//main
