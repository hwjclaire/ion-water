/* The electrostatic energy of NaCl using Ewald sum
 * Copyright 2009-2010 Cheng Zhang
 *
 * This program is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public License can be obtained from 
 * <http://www.gnu.org/licenses/>. */
#include <stdio.h>
#include <math.h>
#define N 2
#define M 7
int main(void)
{
  int i, j, k, nz;
  double x, dg, k0, cwe, ewc = 2.0; /* ewc=sqrt(alpha) */
  double real, recip, self, e;

  real = 0; /* real space sum */
  for (i = 0; i <= N; i++)
  for (j = 0; j <= N; j++)
  for (k = 0; k <= N; k++) {
    nz = (i == 0) + (j == 0) + (k == 0);
    if (nz == 3) continue; /* skip the origin */
    else if (nz == 2) dg = 2;  /* +x, -x */
    else if (nz == 1) dg = 4;  /* +x+y, +x-y, -x+y, -x-y */
    else dg = 8; /* 8 octants */
    if ((i+j+k) % 2 != 0) dg = -dg; /* sign of the charge */
    x = sqrt(i*i + j*j + k*k);
    real += dg * erfc(ewc*x) / x;
  }

  k0 = 2.0*M_PI/2;
  cwe = 0.25*k0*k0/(ewc*ewc);
  recip = 0;
  /* only sum over odd wave vectors, b/c Na cancels Cl */
  for (i = 1; i <= M; i += 2)
  for (j = 1; j <= M; j += 2)
  for (k = 1; k <= M; k += 2) {
    x = i*i + j*j + k*k;
    recip += exp(-x*cwe)/x;
  }
  recip *= 64/(M_PI*2);
  
  self = -2.0*ewc/sqrt(M_PI); /* self energy */
  e = real + recip + self; /* total energy */
  
  printf("real:  %+.14lf\nrecip: %+.14lf\nself:  %+.14lf\nE:     %+.14lf\n",
         real, recip, self, e);
  printf("errors: real: %g, recip: %g\n",
      erfc(ewc*(N+1))/(N+1), exp(-M*M*cwe)/(M*M));
  return 0;
}

