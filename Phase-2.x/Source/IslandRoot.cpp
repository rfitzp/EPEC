// IslandRoot.cpp

#include "Phase.h"

// ############################################################
// Function to find full width in PsiN of magnetic island chain
// ############################################################
double Phase::GetIslandWidth (double A1, double A2, double A3, double Psi)
{
  double Xminus, Xplus;

  GetIslandLimits (A1, A2, A3, Psi, Xminus, Xplus);

  return Xplus - Xminus;
}

// ########################################################################################
// Function to find offset from rational surface in PsiN of center of magnetic island chain
// ########################################################################################
double Phase::GetIslandOffset (double A1, double A2, double A3, double Psi)
{
  double Xminus, Xplus;

  GetIslandLimits (A1, A2, A3, Psi, Xminus, Xplus);

  return (Xplus + Xminus) /2.;
}

// ########################################################
// Function to find limits in PsiN of magnetic island chain 
// ########################################################
void Phase::GetIslandLimits (double A1, double A2, double A3, double Psi, double& Xminus, double& Xplus)
{
  if (fabs (Psi) < 1.e-15)
    {
      Xminus = 0.;
      Xplus  = 0.;

      return;
    }
  
  double LARGE = 1.e6;

  // Calculate basis island width
  double W0 = 4. * sqrt (A1 * fabs (Psi));
  double c = W0*W0/4.;
  
  // Find maximum allowed value of c
  double CMAX = LARGE;
  double cmax = CMAX; 
  if (9.*A2*A2 - 32.*A3 > 0.)
    {
      double xp = (- 3.*A2 + sqrt (9.*A2*A2 - 32.*A3)) /8./A3;
      double xm = (- 3.*A2 - sqrt (9.*A2*A2 - 32.*A3)) /8./A3;

      double fxxp = 2. + 6.*A2*xp + 12.*A3*xp*xp;
      double fxxm = 2. + 6.*A2*xm + 12.*A3*xm*xm;

      double fp, fm;
      if (fxxp < 0.)
	fp = xp*xp + A2*xp*xp*xp + A3*xp*xp*xp*xp;
      else
	fp = CMAX;

      if (fxxm < 0.)
	fm = xm*xm + A2*xm*xm*xm + A3*xm*xm*xm*xm;
      else
	fm = CMAX;

      if (fp < fm)
	cmax = fp;
      else
	cmax = fm;
    }

  // Prevent c from exceeding maximum allowed value
  if (c > cmax)
    c = 0.99*cmax;

  // Solve polynomial equation
  double a[5] = {-c, 0., 1., A2, A3};
  double z[8];

  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc (5);

  gsl_poly_complex_solve (a, 5, w, z);

  gsl_poly_complex_workspace_free (w);

  // Find appropriate real roots of polynomial equation
  int    nroot = 0;
  double x[4];

  for (int i = 0; i < 4; i++)
    if (fabs (z[2*i+1]) < 1.e-15)
      {
	x[i] = z[2*i];
	nroot++;
      }

  if (nroot < 2)
    {
      printf ("FLUX::IslandRoot: Error - nroot = %1d\n", nroot);
      exit (1);
    }

  Xminus = -LARGE;
  Xplus  = +LARGE;
  for (int i = 0; i < nroot; i++)
    if (x[i] < 0.)
      {
	if (x[i] > Xminus)
	  Xminus = x[i];
      }
    else
      {
	if (x[i] < Xplus)
	  Xplus = x[i];
      }
  
  if (Xminus == -LARGE)
    {
      printf ("FLUX::IslandRoot: Error - Xminus = %10.4e\n", Xminus);
      exit (1);
    }
  if (Xplus == LARGE)
    {
      printf ("FLUX::IslandRoot: Error - Xminus = %10.4e\n", Xplus);
      exit (1);
    }
}
