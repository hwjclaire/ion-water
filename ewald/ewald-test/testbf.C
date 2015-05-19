#include <iostream>
#include <vector>
#include "D3vector.h"
#include "Cell.h"

using namespace std;
int main()
{
  const int shell_max = 20;

  double a = 2;

  D3vector atom[8];
  
  double charge[8];
    
  double eps = 0.01;

  std::vector<double> vec;

  vec.resize(9,0);
  vec[0]=2;vec[4]=2;vec[8]=2;

  Cell c(vec);



  atom[0]=D3vector(0,0,0);
  atom[1]=D3vector(1,1,0);
  atom[2]=D3vector(0,1,1);
  atom[3]=D3vector(1,0,1);
  atom[4]=D3vector(1,0,0);
  atom[5]=D3vector(0,1,0);
  atom[6]=D3vector(0,0,1);
  atom[7]=D3vector(1,1,1);

  for ( int i = 0; i < 4; i ++ ) 
  {
    charge[i] = 1;
    charge[i+4] = -1;

  }

  double potential = 0;
  
  int count = 0;
  //loop over layers
  for ( int ishell = 0; ishell < shell_max; ishell ++ )
  {
    for ( int i = - ishell; i <= ishell; i ++ )
    for ( int j = - ishell; j <= ishell; j ++ )
    for ( int k = - ishell; k <= ishell; k ++ )
    {
      if ( fabs(i) < ishell && fabs(j) < ishell && fabs(k) < ishell ) continue;
      //cout << ishell << "    " << i << "   " << j << "   " << k << endl;
      count ++;

      D3vector vec_cell ( i * a, j * a, k * a);



      for ( int iatom = 0; iatom < 8; iatom ++ )
      {
        if ( iatom == 0 && ishell  == 0 ) continue;
    
        D3vector vec1 = atom[iatom] + vec_cell - atom[0];
        //D3vector vec1 = atom[iatom] -atom[0];
        //c.images(vec1);
        //vec1 += vec_cell;

        double dist = length(vec1);
        //if ( dist < eps ) continue;
        //if ( dist > a * ishell  || dist < a * ( ishell -1 ) ) continue;
        potential += charge[iatom] / dist;

        //cout << ishell << "  " << iatom << "   " << dist << endl;

      }

    }
    cout.precision(20);
    cout << ishell << "   " << count << "    "  << potential << endl; 


  }


  cout << count <<endl;














}
