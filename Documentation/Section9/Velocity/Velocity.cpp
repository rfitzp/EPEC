// Program to determine velocity profiles around isolated magnetic island chain

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_ellint.h>

int main ()
{
  FILE* file = fopen("Velocity.out", "w");

  double Xmax = 8.;
  int    N    = 2000;

  for (int i = 0; i <= N; i++)
    {
      double X = - Xmax + (double (i) /double (N)) * 2. * Xmax;
      double V;
      if (fabs(X) <= 2.)
	V = 0.;
      else
	V = (M_PI/2.) /gsl_sf_ellint_Ecomp (2./fabs(X), GSL_PREC_DOUBLE);

      fprintf (file, "%e %e\n", X, V);
    }
  
  fclose (file);
}
