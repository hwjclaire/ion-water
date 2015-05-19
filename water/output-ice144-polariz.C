// output-ice.x
// For a 144-D2O ice Ih surface simulation ( 3 bilayers )
// and compute molecular polarizabilities using ewald sum
// Quan (Andy) Wan Sat Oct  5 00:39:22 PDT 2013

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
    cout << "output-ice.x [input xyz] [input mlwf] [input quadrupole] [input polar] #frame nMO nskip\n";
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
  foxyz.open("out.xyz");

  Stat quad(-5,5,0.01);
  Stat dip(0,2,0.01);
  Stat owf(0,2,0.01);

  ofstream fdip,fq,fowf;
  fdip.open("statdip");
  fq.open("statquad");
  fowf.open("statowf");

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
        name << "1water" << i + 1 << ".dat";
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

      // output xyz 
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
        owf.add(min_dist);
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
        for ( int k = 0; k < 3 ; k ++ )
          quad.add ( waterset[i].quad_eigval()[k] );
        dip.add(length(waterset[i].dipole()));
  
        //output each molecules
        //polarmlwf[i] << waterset[i].dipole();
        //waterset[i].polariz().print_inline(polarmlwf[i]);
        //polarmlwf[i] << endl;
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

      tm_alpha.start();
      polariz -> Calculate_mol_Ewald();
      tm_alpha.stop();


      tm_output.start();

#if 1
      D3vector d1(0,0,0),d2(0,0,0),d3(0,0,0),d4(0,0,0),d5(0,0,0),d6(0,0,0);
      Tensor p1,p2,p3,p4,p5,p6;
      int count1=0,count2=0,count3=0,count4=0,count5=0,count6=0;

      for ( int iwater = 0; iwater < waterset.size(); iwater ++ )
      {


        //Tensor p(&(polariz->mol_elocal()[iwater][0]));
        Tensor p(waterset[iwater].alpha_tensor());
        polarmlwf[iwater] << waterset[iwater].dipole();
        p.print_inline(polarmlwf[iwater]);
        polarmlwf[iwater] << endl;


        double z = waterset[iwater].o()->x().z;
        //cout << z << endl;
        if ( z < 6 ) 
        {
          d1 += waterset[iwater].dipole();
          p1 += Tensor(&(polariz->mol_elocal()[iwater][0]));
          count1++;
        }
        else if ( z < 15 )
        {
          d2 += waterset[iwater].dipole();
          p2 += Tensor(&(polariz->mol_elocal()[iwater][0]));
          count2++;
        }
        else if ( z < 21 )
        {
          d3 += waterset[iwater].dipole();
          p3 += Tensor(&(polariz->mol_elocal()[iwater][0]));
          count3++;
        }
        else if ( z < 27 ) 
        {
          d4 += waterset[iwater].dipole();
          p4 += Tensor(&(polariz->mol_elocal()[iwater][0]));
          count4++;
        }
        else if ( z < 34 )
        {
          d5 += waterset[iwater].dipole();
          p5 += Tensor(&(polariz->mol_elocal()[iwater][0]));
          count5++;
        }
        else
        {
          d6 += waterset[iwater].dipole();
          p6 += Tensor(&(polariz->mol_elocal()[iwater][0]));
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


      assert( count1 == 24 );
      assert( count2 == 24 );
      assert( count3 == 24 );
      assert( count4 == 24 );
      assert( count5 == 24 );
      assert( count6 == 24 );

      fbo << d1;
      p1.print_inline(fbo);
      fbo << endl;

      fbo2 << d2;
      p2.print_inline(fbo2);
      fbo2 << endl;

      fbo3 << d3;
      p3.print_inline(fbo3);
      fbo3 << endl;

      fup << d6;
      p6.print_inline(fup);
      fup << endl;

      fup2 << d5;
      p5.print_inline(fup2);
      fup2 << endl;

      fup3 << d4;
      p4.print_inline(fup3);
      fup3 << endl;


      Tensor p = p3+p4;
      fmid << d3+d4;
      p.print_inline(fmid);
      fmid << endl;



      fdipole << d1+d2+d3+d4+d5+d6 << endl;

#endif
      tm_output.stop();

    }// if iframe % nskip

  }// for iframe

  quad.print(fq);
  dip.print(fdip);
  owf.print(fowf);

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
