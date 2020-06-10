// ToroidalP.cpp

#include "Flux.h"

// ################################################################
// Function to return associated Legendre function
//  
//   P^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// Routine sums hypergeometric series representation of function to 
// accuracy of about 1 in 10^12.
//
// Reference: Bateman Manuscript Project, Vol. I, p. 124, Eq. (16)
// ################################################################

double Flux::ToroidalP (int m, int n, double z)
{
  // Check argument
  if (z < 1.)
    {
      printf ("FLUX:: ToroidalP: Error z < 1.\n");
      exit (1);
    }	

  // Calculate factor multiplying hypergeometric series
  double x  = (z - 1.) /(z + 1.);
  double d  = pow (0.5 + z/2., double (abs (n)) - 0.5);
  d        *= pow (x, double (abs (m)) /2.);
  
  if (m > 0)
    for (int j = 1; j <= m; j++)
      d *= (double (n*n - j*j + j) - 0.25) / double (j);
  else if (m < 0)
    {
      int mp = - m;
      for (int j = 1; j <= mp; j++)
	d /= double (j);
    }

  // Sum hypergeometric series
  double a, c, b, e, f, q, r;
  a = 0.5 - double (abs (n));
  c = 1.  + double (abs (m));
  b = a + c - 1.;
  e = 1.;
  f = 1.;
  q = 0.;
  
  do 
    {
      e *= x * (a + q) * (b + q) /(c + q) /(1. + q);
      f += e;
      q += 1.;
      
      r  = (a + q) /(a + q - 1.);
      r *= (b + q) /(b + q - 1.);
      r *= (c + q - 1.) /(c + q);
      r *= q /(q + 1.);
    }
  while (fabs (e) > 1.e-15 || fabs (r) > 1./x || q < - double (n));
  
  return f * d;
} 
  
