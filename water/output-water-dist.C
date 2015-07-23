// output-water-mlwf.x
// For simulation of bulk water/ice
// output mlwf spread and polarizability
// Quan (Andy) Wan Sep 26 12:07:11 CDT 2014

#include <iostream>
#include <iomanip>
#include <vector>
//#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
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
    cout << "output-water-mlwf.x [input xyz] [input mlwf] [input quadrupole] [input polar] #frame nwater nskip\n";
    return 1;
  }
 
  ifstream fxyz, fmlwf, fquad,fpolar;
  fxyz.open(argv[1]);
  fmlwf.open(argv[2]);
  fquad.open(argv[3]);
  fpolar.open(argv[4]);

  typedef std::map<double,Water*>  DistMap;



  int nframe = atoi(argv[5]), nwater = atoi(argv[6]), nskip = atoi(argv[7]), natom;
  int nmo = nwater * 4;
  double volume;
  Cell c;
  std::vector<Atom> atomset;
  std::vector<Mlwf> mlwfset;
  std::vector<Water> waterset;
  std::vector<Mol*> molset;

  //ofstream polarmlwf[nwater];


  Stat quad(-5,5,0.01);
  Stat dip(0,2,0.01);

  ofstream fdip,fq;
  //fdip.open("statdip");
  //fq.open("statquad");

  ofstream ofmlwf, ofwater;
  //ofmlwf.open("output-mlwf");
  ofwater.open("output-water");


  Polariz *polariz;

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
        name << "waterquad" << i + 1 << ".dat";
        string filename = name.str();
        char * fptr = (char*)filename.c_str();
        //polarmlwf[i].open(fptr);
        //polarmlwf[i].setf(ios::fixed, ios::floatfield);
        //polarmlwf[i].setf(ios::right, ios::adjustfield);
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

      // compute mlwf every 500 step
      if ( iframe % 100 != 0 ) continue;
      //cout << iframe << endl;
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
        //cout << min->number() << endl;
        min -> add_wf(*pwf);
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
        waterset[i].getaxis();
        //waterset[i].Compute_quad();
        waterset[i].Compute_polariz();
        //for ( int k = 0; k < 3 ; k ++ )
        //  quad.add ( waterset[i].quad_eigval()[k] );
        //dip.add(length(waterset[i].dipole()));

      }

      if ( iframe == 0 )
      {
        molset.resize(waterset.size());
        for ( int i = 0; i < molset.size(); i ++ )
          molset[i] = (Mol*)&waterset[i];

        polariz = new Polariz(c, molset, mlwfset);
      }

      tm_alpha.start();
      polariz -> Calculate_mol_Ewald();
      //polariz -> Calculate_mol();
      tm_alpha.stop();

      //store oxygen oxygen distance between water
      std::vector<double> dist(nwater*nwater);
  
      // compute distance matrix and store in dist
      for ( int i = 0; i < waterset.size(); i ++ ) 
      {

        //cout << i << endl;
        dist[i*nwater+i] = 0;

        waterset[i].rotate_alpha_tensor();
        waterset[i].rotate_polar_tensor();

        for ( int j = 0; j < waterset.size(); j ++ )
        {
          if ( i == j ) continue;
          double d0 = waterset[j].o()->distance(waterset[i].o(),&c);

          dist[i*nwater+j]=d0;
          dist[j*nwater+i]=d0;
        }

        //ofwater << waterset[i].polariz().trace()/3 << " " ;
        //ofwater << Tensor(molset[i]->alpha_tensor()).trace()/3 << endl;


      }

      for ( int i = 0; i < waterset.size(); i ++ )
      {

        // build distance map to water pointer
        // this is used to sort water from near to far
        DistMap dmap;
        for ( int j = 0; j < waterset.size(); j ++ )
          dmap[dist[i*nwater+j]] = &waterset[j];

        // output inner distances
        D3vector r01 = waterset[i].oh(0);
        D3vector r02 = waterset[i].oh(1);
        double sr = length(r01)+length(r02);
        double dr = length(r01)-length(r02);
        ofwater << sr << " " << dr*dr << " " << r01*r02 << " ";

        Water* iwater =  &waterset[i];
        // output distances with neighboring waters
        for ( DistMap::iterator j = dmap.begin(); j != dmap.end(); j ++ )
        {
          Water* jwater =  j->second;
          if ( jwater == iwater ) continue;   
          //output inner diatances

          //output distances between two waters
          // OO distance
          ofwater << (jwater->o()->distance(iwater->o(),&c)) << " ";
          // OH distance
          double r1 = (jwater->o()->distance(iwater->h(0),&c));
          double r2 = (jwater->o()->distance(iwater->h(1),&c));
          ofwater << ( r1 < r2 ? r1 : r2 ) << " " << ( r1 < r2 ? r2 : r1 ) << " ";

          // HO distance
          r1 = (jwater->h(0)->distance(iwater->o(),&c));
          r2 = (jwater->h(1)->distance(iwater->o(),&c));

          // HH distances
          if ( r1 < r2 )
          {
            double r3 = (jwater->h(0)->distance(iwater->h(0),&c));
            double r4 = (jwater->h(0)->distance(iwater->h(1),&c));
            ofwater << r1 << " " << ( r3 < r4 ? r3 : r4 ) << " " << ( r3 < r4 ? r4 : r3 ) << " ";

            r3 = (jwater->h(1)->distance(iwater->h(0),&c));
            r4 = (jwater->h(1)->distance(iwater->h(1),&c));
            ofwater << r2 << " " << ( r3 < r4 ? r3 : r4 ) << " " << ( r3 < r4 ? r4 : r3 ) << " ";

          }
          else
          {
            double r3 = (jwater->h(1)->distance(iwater->h(0),&c));
            double r4 = (jwater->h(1)->distance(iwater->h(1),&c));
            ofwater << r2 << " " << ( r3 < r4 ? r3 : r4 ) << " " << ( r3 < r4 ? r4 : r3 ) << " ";

            r3 = (jwater->h(0)->distance(iwater->h(0),&c));
            r4 = (jwater->h(0)->distance(iwater->h(1),&c));
            ofwater << r1 << " " << ( r3 < r4 ? r3 : r4 ) << " " << ( r3 < r4 ? r4 : r3 ) << " ";

          }


        } // for j

        ofwater << length(iwater->dipole()) << " ";
        ofwater << iwater->dipole() * D3vector(iwater->axis_ptr(0)) << " ";
        ofwater << iwater->dipole() * D3vector(iwater->axis_ptr(1)) << " ";
        ofwater << iwater->dipole() * D3vector(iwater->axis_ptr(2)) << " ";
        ofwater << iwater->alpha_tensor_axis(0) << " ";
        ofwater << iwater->alpha_tensor_axis(4) << " ";
        ofwater << iwater->alpha_tensor_axis(8) << " ";
        ofwater << iwater->alpha_tensor_axis(1) << " ";
        ofwater << iwater->alpha_tensor_axis(3) << " ";
        ofwater << iwater->alpha_tensor_axis(2) << " ";
  
        ofwater << iwater->polar_tensor_axis(0) << " ";
        ofwater << iwater->polar_tensor_axis(4) << " ";
        ofwater << iwater->polar_tensor_axis(8) << " ";
        ofwater << iwater->polar_tensor_axis(1) << " ";
        ofwater << iwater->polar_tensor_axis(3) << " ";
        ofwater << iwater->polar_tensor_axis(2) << " ";
        ofwater << endl;
             


      } // for i


      tm_misc.stop();
      //cout << "Check Done" << endl;
  
      tm_output.start();
      for ( int i = 0; i < waterset.size(); i ++ )
      {
        
        //polarmlwf[i] << waterset[i].dipole();
        //waterset[i].quad().print_inline(polarmlwf[i]);
        //polarmlwf[i] << endl;

      }

      tm_output.stop();

    }// if iframe % nskip

  }// for iframe

  quad.print(fq);
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
  fquad.close();

}//main
