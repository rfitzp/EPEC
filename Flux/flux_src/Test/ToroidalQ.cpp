// ToroidalQ.cpp

#include "Flux.h"

// ################################################################
// Function to return associated Legendre function
//  
//   Q^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// Routine sums hypergeometric series representation of function to 
// accuracy of about 1 in 10^12.
//
// Reference: Bateman Manuscript Project, Vol. I, p. 134, Eq. (41)
// ################################################################
double Flux::ToroidalQ (int m, int n, double z)
{
  // Check argument
  if (z < 1.)
    {
      printf ("Error: ToroidalQ: z < 1.\n");
      exit (1);
    }

  // Calculate factor multiplying hypergeometric series
  double y  = sqrt ((z*z - 1.) /z/z);
  double d  = M_PI / sqrt (2.*z); 
  int    na = abs (n);
  if (m < 0)
    {
      int mp = - m;
      for (int j = 1; j <= mp; j++)
	d /= double (n*n - j*j + j) - 0.25;
    }
  int ma = abs (m);
  if (na > 0)
    for (int k = 1; k <= na; k++)
      d *= (double (ma + k) - 0.5) /2./z / double (k);
  if (ma > 0)
    for (int l = 1; l <= ma; l++)
      d *= - y * (double (l) - 0.5);

  // Sum hypergeometric series
  double a, b, c, e, f, q, r;
  a = 0.5 * (1.5 + double (ma + na));
  b = a - 0.5;
  c = double (na + 1);
  e = 1.;
  f = 1.;
  q = 0.;
  
  do
    {
      e *= (a + q) * (b + q) / (c + q) / (1. + q) / z / z;
      f += e;
      q += 1.;
      
      r  = (a + q) / (a + q - 1.);
      r *= (b + q) / (b + q - 1.);
      r *= (c + q - 1.) / (c + q);
      r *= q / (q + 1.);
    }
  while (e > 1.e-15 || r > z*z);
  
  return f * d;
} 
  
