// output-ice.x
// For a 96-D2O and alumina interface simulation
// Quan (Andy) Wan Thu Feb  6 15:55:27 PST 2014

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
  std::vector<Water> alohset1, alohset2;
  std::vector<Mol*> molset;

  ofstream alo, alo1, alo2;
  alo.open("alo.dat");
  alo1.open("alos1.dat");
  alo2.open("alos2.dat");

  ofstream waterox[96];
  ofstream polarmlwf[96];
  ofstream polarmlwf1[36];
  ofstream polarmlwf2[36];
  ofstream watermlwf,watermlwfpolar;
  watermlwf.open("watermlwf.dat");
  watermlwfpolar.open("watermlwfpolar.dat");

  // list of water within x a.u. of the alumina surface.
  // surface 2 is on top of surface 1 if look in the z direction
  ofstream wlist1[18], wlist2[18];
  
  for ( int i = 0; i < 18; i ++ )
  {
    ostringstream name;
    name << "s1-" << i + 1 << "au.dat";
    string filename = name.str();
    char * fptr = (char*)filename.c_str();
    wlist1[i].open(fptr);
  }
  for ( int i = 0; i < 18; i ++ )
  {
    ostringstream name;
    name << "s2-" << i + 1 << "au.dat";
    string filename = name.str();
    char * fptr = (char*)filename.c_str();
    wlist2[i].open(fptr);
  }

  Stat quad(-5,5,0.01);
  Stat dip(0,2,0.01);
  Stat dist_alo(0,28,0.02);
  Stat angle_aloh(-2,2,0.01);

  std::vector < double > angle_water_sum;
  std::vector < int > angle_water_count;

  ofstream fdip,fq,fdist,fanglealoh,fanglewater;
  fdip.open("statdip");
  fq.open("statquad");
  fdist.open("statdist");
  fanglealoh.open("statanglealoh");
  fanglewater.open("statanglewater");

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
    tm_read.start();
    for ( int iatom = 0; iatom < natom; iatom ++ )
    {
      string name;

      fxyz >> name;

      if ( iframe == 0 )
      {
        if ( name.compare("OW") == 0 || name.compare("Oh") == 0 || name.compare("Ox") == 0 )
        {
          
          atomset[iatom].setmass(15.999);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(6.0);
          
        }
        else if ( name.compare("HW") == 0 || name.compare("Hh") == 0 )
        {
          atomset[iatom].setmass(2.014);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(1.0);
        }
        else if ( name.compare("Al") == 0 )
        {
          atomset[iatom].setmass(26.981);
          atomset[iatom].setname(name);
          atomset[iatom].setcharge(1.0);
        }
        else
        {
          return 1;
        }
      }//if iframe == 0
      // read position, velocity
      fxyz >> atomset[iatom];
      //cout << atomset[iatom].x()<< endl;

    }//for iatom
    tm_read.stop();

    //computing the dividing surface.
    double z0=0;
    double alozmin = z0;
    double alozmax = z0;
    if ( nskip == 0 || iframe % nskip == 0 )
    {
      z0=0;
      int n0=0;
      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        if ( atomset[iatom].name().compare("Ox") == 0 )
        {
          z0 += atomset[iatom].x().z;
          n0++;
        }
      }
      z0 /= n0;
  
      alozmin = z0;
      alozmax = z0;
  
      for ( int iatom = 0; iatom < natom; iatom ++ ) 
      { 
        if ( atomset[iatom].name().compare("Hh") == 0 ) 
        { 
          double z = atomset[iatom].x().z;
          if ( z < alozmin ) alozmin = z;
          if ( z > alozmax ) alozmax = z;
        } 
      } 
  
      alo << z0 << " " << alozmin << " " << alozmax << endl;
    } // if nskip

    
    //assign H to O
    if ( iframe == 0 )  //assign molecules 
    {

      int watercount = 0;
      int alohcount1 = 0;
      int alohcount2 = 0;
 
      for ( int i = 0; i < 96; i ++ )
        waterset.push_back(Water(i,c));
      for ( int i = 0; i < 36; i ++ )
        alohset1.push_back(Water(i,c));
      for ( int i = 0; i < 36; i ++ )
        alohset2.push_back(Water(i,c));

 
      for ( int iatom = 0; iatom < natom; iatom ++ )
      {
        if ( atomset[iatom].name().compare("OW") == 0 )
        {
          waterset[watercount].add_oxygen ( atomset[iatom] );
          watercount ++;
        }
        else if ( atomset[iatom].name().compare("Oh") == 0 )
        {
          if ( atomset[iatom].x().z < z0 )
          {
            alohset1[alohcount1++].add_oxygen ( atomset[iatom] );
          }
          else
          {
            alohset2[alohcount2++].add_oxygen ( atomset[iatom] );
          }

        }
      }

      cout << "watercount=" << watercount << endl;
      cout << "alohcount1=" << alohcount1 << endl;
      cout << "alohcount2=" << alohcount2 << endl;

      //assign water hydrogen to water oxygen
      for ( int iatom = 0; iatom < natom; iatom ++)
      { 
        Atom * patom = & atomset[iatom];
        if ( patom -> name().compare("HW") ) continue;//is not H
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

      //assert ( waterset.size() * 4 == nmo );
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        assert(waterset[i].isoh2());
        waterset[i].check_atom();
      }
      cout << "done add water h" << endl;


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

      for ( int i = 0; i < waterset.size(); i ++ )
      {
        ostringstream name;
        name << "waterox" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        waterox[i].open(fptr);
        waterox[i].setf(ios::fixed, ios::floatfield);
        waterox[i].setf(ios::right, ios::adjustfield);
        waterox[i].precision(4);
    
      }

      for ( int i = 0; i < alohset1.size(); i ++ )
      {
        ostringstream name;
        name << "aloh1" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmlwf1[i].open(fptr);
        polarmlwf1[i].setf(ios::fixed, ios::floatfield);
        polarmlwf1[i].setf(ios::right, ios::adjustfield);
        polarmlwf1[i].precision(4);
      }

      for ( int i = 0; i < alohset2.size(); i ++ )
      {
        ostringstream name;
        name << "aloh2" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        polarmlwf2[i].open(fptr);
        polarmlwf2[i].setf(ios::fixed, ios::floatfield);
        polarmlwf2[i].setf(ios::right, ios::adjustfield);
        polarmlwf2[i].precision(4);
      }




      cout << "done add wf to water" << endl;

      //assign Al-OH hydrogen to Al- oxygen
      for ( int iatom = 0; iatom < natom; iatom ++)
      {
        Atom * patom = & atomset[iatom];
        if ( patom -> name().compare("Hh") ) continue;//is not H
        double min_dist = 1000; Mol * min;
        //look in surface1
        for ( int iwater = 0; iwater < alohset1.size(); iwater ++)
        {
          Mol* pwater = &alohset1[iwater];
          double dist = pwater -> distance(patom);
          if ( dist < min_dist)
          {
            min = pwater;
            min_dist = dist;
          }//if

        }//for jatom
        //look in surface2
        for ( int iwater = 0; iwater < alohset2.size(); iwater ++)
        {
          Mol* pwater = &alohset2[iwater];
          double dist = pwater -> distance(patom);
          if ( dist < min_dist)
          {
            min = pwater;
            min_dist = dist;
          }//if

        }//for jatom

        min -> add_hydrogen ( *patom );
      }//for iatom, done adding hydrogen to Al-oxygen
      
      for ( int i = 0; i < alohset1.size(); i ++ )
      {
        alohset1[i].check_atom();
        alohset2[i].check_atom();
        assert(alohset1[i].isoh());
        assert(alohset2[i].isoh());
      }
      cout << "addh\n";

    }// if ifram == 0;  

    if ( nskip == 0 || iframe % nskip == 0 )
    {

      // output list of water close to surface
      // output water dipole angle distribution
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        double z = waterset[i].o()->x().z;
        double zc = c.z();
        double dist1 = ( z - alozmax ) / zc;
        dist1 = ( dist1 - floor(dist1) ) * zc;
  
        double dist2 = ( alozmin - z ) / zc;
        dist2 = ( dist2 - floor(dist2) ) * zc;
  
        //cout << i << " " << dist1 << " " << dist2 << endl;

        dist_alo.add(dist1);

        for ( int j = 0; j < 18; j ++ )
        {
          if ( dist1 < (1+j) )
            wlist1[j] << i << " ";
          if ( dist2 < (1+j) )
            wlist2[j] << i << " ";
        }
  
      }
      for ( int j = 0; j < 18; j ++ )
      {
        wlist1[j] << endl;
        wlist2[j] << endl;
      }


      //add mlwfs
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
      for ( int i = 0; i < alohset1.size(); i ++ )
      {
        //update oh distances.
        alohset1[i].watermove();
        alohset1[i].reset_wf();
      }
      for ( int i = 0; i < alohset2.size(); i ++ )
      {
        //update oh distances.
        alohset2[i].watermove();
        alohset2[i].reset_wf();
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
          //look in waters
          Water * pwater = &waterset[iwater];
          //if (waterset[iwater].wf_full()) continue;
          double dist = pwater -> distance(pwf);
          if ( dist < min_dist )
          {
            min = pwater;
            min_dist = dist;
          }
        }
        for ( int iwater = 0; iwater < alohset1.size(); iwater ++)
        {
          //look in alumina surface OH
          Water * pwater = &alohset1[iwater];
          //if (waterset[iwater].wf_full()) continue;
          double dist = pwater -> distance(pwf);
          if ( dist < min_dist )
          {
            min = pwater;
            min_dist = dist;
          }
        }
        for ( int iwater = 0; iwater < alohset2.size(); iwater ++)
        {
          //look in alumina surface OH
          Water * pwater = &alohset2[iwater];
          //if (waterset[iwater].wf_full()) continue;
          double dist = pwater -> distance(pwf);
          if ( dist < min_dist )
          {
            min = pwater;
            min_dist = dist;
          }
        }

        if ( min_dist > 1.5 ) continue;
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
        //waterset[i].Compute_polariz();
        //for ( int k = 0; k < 3 ; k ++ )
        //  quad.add ( waterset[i].quad_eigval()[k] );
        dip.add(length(waterset[i].dipole()));
      }

      for ( int i = 0; i < alohset1.size(); i ++ )
      {
        alohset1[i].mlwfmove();
        alohset1[i].check_wf();
        alohset1[i].Compute_cc();
        alohset1[i].Compute_cm();
        alohset1[i].Compute_dipole();
        alohset1[i].Compute_polariz();
        angle_aloh.add(angle(alohset1[i].oh(0),D3vector(0,0,1)));
      }

      for ( int i = 0; i < alohset2.size(); i ++ )
      {
        alohset2[i].mlwfmove();
        alohset2[i].check_wf();
        alohset2[i].Compute_cc();
        alohset2[i].Compute_cm();
        alohset2[i].Compute_dipole();
        alohset2[i].Compute_polariz();
        angle_aloh.add(angle(alohset2[i].oh(0),D3vector(0,0,-1)));
      }


      angle_water_sum.resize(400,0);
      angle_water_count.resize(400,0);
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        double z = waterset[i].o()->x().z;
        double zc = c.z();
        double dist1 = ( z - alozmax ) / zc;
        dist1 = ( dist1 - floor(dist1) ) * zc;

        double dist2 = ( alozmin - z ) / zc;
        dist2 = ( dist2 - floor(dist2) ) * zc;

        for ( int j = 0; j < 400; j ++ )
        {
          if ( dist1 < zc / 400 * j )
          {
            //cout << waterset[i].dipole() << "  "
            //     << angle ( waterset[i].dipole(), D3vector(0,0,1) ) << endl;
            angle_water_sum[j] += angle (
                          waterset[i].dipole(), D3vector(0,0,1) );
            angle_water_count[j] ++;
            break;
          }
        }

      }


      tm_misc.stop();
      //cout << "Check Done" << endl;
  
      tm_output.start();


      //output dipole, polarizability for each individual water and OH groups
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        polarmlwf[i] << waterset[i].dipole();
        waterset[i].polariz().print_inline(polarmlwf[i]);
        polarmlwf[i] << endl;
        waterox[i] << waterset[i].o()->x() << endl;

        for ( int j = 0; j < 4; j ++ )
        {
          watermlwf << waterset[i].wf(j)->x() << " " << waterset[i].wf(j)->spread() << endl;
          for ( int k = 0; k < 6; k ++ )
            watermlwfpolar << waterset[i].wf(j)->polar()[k] << " ";
          watermlwfpolar << endl;
        }
      }

      D3vector d1(0,0,0), d2(0,0,0);
      Tensor p1, p2;

      for ( int i = 0; i < alohset1.size(); i ++ )
      {
        polarmlwf1[i] << alohset1[i].dipole();
        alohset1[i].polariz().print_inline(polarmlwf1[i]);
        polarmlwf1[i] << endl;
        d1 += alohset1[i].dipole();
        p1 += alohset1[i].polariz();
      }

      alo1 << d1;
      p1.print_inline(alo1);
      alo1 << endl;

      for ( int i = 0; i < alohset2.size(); i ++ )
      {
        polarmlwf2[i] << alohset2[i].dipole();
        alohset2[i].polariz().print_inline(polarmlwf2[i]);
        polarmlwf2[i] << endl;
        d2 += alohset2[i].dipole();
        p2 += alohset2[i].polariz();
      }

      alo2 << d2;
      p2.print_inline(alo2);
      alo2 << endl;


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

  for ( int i = 0; i < waterset.size(); i ++ ) polarmlwf[i].close();
  for ( int i = 0; i < alohset1.size(); i ++ ) polarmlwf1[i].close();
  for ( int i = 0; i < alohset2.size(); i ++ ) polarmlwf2[i].close();

  //quad.print(fq);
  dip.print(fdip);
  dist_alo.print(fdist);

  for ( int i = 0; i < 400; i ++ )  
  {
    if ( angle_water_count[i] == 0 )
      fanglewater << i*c.z()/400 << "  " << 0 << "  " << 0 << endl;
    else
      fanglewater << i*c.z()/400 << "  " << angle_water_count[i]
                  << "  " << angle_water_sum[i]/angle_water_count[i] << endl;
  }


  fxyz.close();
  fmlwf.close();
  fpolar.close();
  //fquad.close();

}//main
