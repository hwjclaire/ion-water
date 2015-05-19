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
  fdip.open("statdip");
  fq.open("statquad");

  ofstream ofmlwf, ofwater;
  ofmlwf.open("output-mlwf");
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
      if ( iframe % 500 != 0 ) continue;

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
        min -> add_wf(*pwf,false);
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

      for ( int i = 0; i < waterset.size(); i ++ ) 
      {


        if ( iframe % 500 == 0 )
        {

          // find closest water
          Water * nn0, * nn1, *nn2, *nn3;
          double dist0 = 100, dist1 = 100, dist2=100, dist3=100;
          for ( int j = 0; j < waterset.size(); j ++ )
          {
            if ( i == j ) continue;
            double d0 = waterset[j].o()->distance(waterset[i].wf(0),&c);
            double d1 = waterset[j].o()->distance(waterset[i].wf(1),&c);
            double d2 = waterset[j].o()->distance(waterset[i].wf(2),&c);
            double d3 = waterset[j].o()->distance(waterset[i].wf(3),&c);
            if ( d0 < dist0 )
            {
              nn0 = &waterset[j];
              dist0 = d0;
            }
  
            if ( d1 < dist1 )
            {
              nn1 = &waterset[j];
              dist1 = d1;
            }

            if ( d2 < dist2 )
            {
              nn2 = &waterset[j];
              dist2 = d2;
            }

            if ( d3 < dist3 )
            {
              nn3 = &waterset[j];
              dist3 = d3;
            }



          }
          //cout << i << " " << nn0->number() << " " << nn1-> number() << " ";
          //cout << i << " " << nn2->number() << " " << nn3-> number() << endl;
          //cout << i << " " << dist0 << " " << dist1 << " ";
          //cout << i << " " << dist2 << " " << dist3 << endl;


          ofmlwf << waterset[i].wf(0)->spread() << " " << waterset[i].wf(1)->spread() << " ";
          ofmlwf << waterset[i].wf(2)->spread() + waterset[i].wf(3)->spread() << " ";
          ofmlwf << length(waterset[i].owf(0)) << " " << length(waterset[i].owf(1)) << " ";
          ofmlwf << length(waterset[i].owf(2)) + length(waterset[i].owf(3)) << " ";
          ofmlwf << length(waterset[i].oh(0)) << " " << length(waterset[i].oh(1)) << " ";
          ofmlwf << angle(waterset[i].oh(0),waterset[i].oh(1)) << " ";
          ofmlwf << angle(waterset[i].oh(0),waterset[i].owf(0)) << " ";
          // nearest water
          ofmlwf << dist0 << " " << dist1 << " " << dist2 << " " << dist3 <<  " ";
          ofmlwf << nn0->o()->distance(waterset[i].h(0),&c) << " ";
          ofmlwf << angle(nn0->o()->x()-waterset[i].h(0)->x(),waterset[i].oh(0)) << " ";
          ofmlwf << nn1->o()->distance(waterset[i].h(1),&c) << " ";
          ofmlwf << angle(nn1->o()->x()-waterset[i].h(1)->x(),waterset[i].oh(1)) << " ";


          ofmlwf << waterset[i].wf(0)->polariz().trace()/3 << endl;

          ofmlwf << waterset[i].wf(1)->spread() << " " << waterset[i].wf(0)->spread() << " ";
          ofmlwf << waterset[i].wf(2)->spread() + waterset[i].wf(3)->spread() << " ";
          ofmlwf << length(waterset[i].owf(1)) << " " << length(waterset[i].owf(0)) << " ";
          ofmlwf << length(waterset[i].owf(2)) + length(waterset[i].owf(3)) << " ";
          ofmlwf << length(waterset[i].oh(1)) << " " << length(waterset[i].oh(0)) << " ";
          ofmlwf << angle(waterset[i].oh(0),waterset[i].oh(1)) << " ";
          ofmlwf << angle(waterset[i].oh(1),waterset[i].owf(1)) << " ";
          // nearest water
          ofmlwf << dist1 << " " << dist0 << " " << dist2 << " " << dist3 <<  " ";
          ofmlwf << nn1->o()->distance(waterset[i].h(1),&c) << " ";
          ofmlwf << angle(nn1->o()->x()-waterset[i].h(1)->x(),waterset[i].oh(1)) << " ";
          ofmlwf << nn0->o()->distance(waterset[i].h(0),&c) << " ";
          ofmlwf << angle(nn0->o()->x()-waterset[i].h(0)->x(),waterset[i].oh(0)) << " ";
          ofmlwf << waterset[i].wf(1)->polariz().trace()/3 << endl;

          ofwater << waterset[i].wf(0)->spread() + waterset[i].wf(1)->spread() +
                     waterset[i].wf(2)->spread() + waterset[i].wf(3)->spread() << " ";
          ofwater << length(waterset[i].owf(0)) << " " << length(waterset[i].owf(1)) << " ";
          ofwater << length(waterset[i].owf(2)) << " " << length(waterset[i].owf(3)) << " ";
          ofwater << angle(waterset[i].oh(0),waterset[i].owf(0)) << " ";
          ofwater << angle(waterset[i].oh(1),waterset[i].owf(1)) << " ";
          ofwater << length(waterset[i].oh(0)) << " " << length(waterset[i].oh(1)) << " ";
          ofwater << angle(waterset[i].oh(0),waterset[i].oh(1)) << " ";

          // 1st neighbour
          ofwater << dist0 << " " << dist1 << " " << dist2 << " " << dist3 <<  " ";
          ofwater << nn0->o()->distance(waterset[i].h(0),&c) << " ";
          ofwater << angle(nn0->o()->x()-waterset[i].h(0)->x(),waterset[i].oh(0)) << " ";
          ofwater << length(nn0->oh(0)) << " " << length(nn0->oh(1)) << " ";
          ofwater << angle(nn0->oh(0),nn0->oh(1)) << " ";

          ofwater << nn1->o()->distance(waterset[i].h(1),&c) << " ";
          ofwater << angle(nn1->o()->x()-waterset[i].h(1)->x(),waterset[i].oh(1)) << " ";
          ofwater << length(nn1->oh(0)) << " " << length(nn1->oh(1)) << " ";
          ofwater << angle(nn1->oh(0),nn1->oh(1)) << " ";

          ofwater << nn2->o()->distance(waterset[i].o(),&c) << " ";
          ofwater << angle(nn0->o()->x()-waterset[i].o()->x(),waterset[i].owf(2)) << " ";
          ofwater << length(nn2->oh(0)) << " " << length(nn2->oh(1)) << " ";
          ofwater << angle(nn2->oh(0),nn2->oh(1)) << " ";

          ofwater << nn3->o()->distance(waterset[i].o(),&c) << " ";
          ofwater << angle(nn1->o()->x()-waterset[i].o()->x(),waterset[i].owf(3)) << " ";
          ofwater << length(nn3->oh(0)) << " " << length(nn3->oh(1)) << " ";
          ofwater << angle(nn3->oh(0),nn3->oh(1)) << " ";
          


          ofwater << waterset[i].polariz().trace()/3 << " " ;
          ofwater << Tensor(molset[i]->alpha_tensor()).trace()/3 << endl;


          




        }

      }

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
