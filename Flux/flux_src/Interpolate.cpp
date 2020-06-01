// Interpolate.cpp

#include "Flux.h"

// ##############################################
// 1D interpolation function with nonuniform grid
// order = 0: Y(x)
// order = 1: dY/dx
// ##############################################
double Flux::Interpolate (int I, double* X, double* Y, double x, int order)
{  
  int i0 = 0;
  for (int i = 1; i < I; i++)
    if (x > X[i])
      i0 = i;
  if (i0 < 0 || i0 > I-1)
    {
      printf ("FLUX::Interpolate - Interpolation error: I = %3d i0 = %3d x = %11.4e X[0] = %11.4e X[I-1] = %11.4e\n",
	      I, i0, x, X[0], X[I-1]);
      exit (1);
    }
  if (x - X[i0] > 0.5 * (X[i0+1] - X[i0]))
    i0 += 1;
  if (i0 == 0)
    i0 += 1;
  if (i0 == I-1)
    i0 -= 1;

  double val;
  if (order == 0)
    {
      double sm = (x - X[i0  ]) * (x - X[i0+1]) /(X[i0-1] - X[i0  ]) /(X[i0-1] - X[i0+1]);
      double s0 = (x - X[i0-1]) * (x - X[i0+1]) /(X[i0  ] - X[i0-1]) /(X[i0  ] - X[i0+1]);
      double s1 = (x - X[i0-1]) * (x - X[i0  ]) /(X[i0+1] - X[i0-1]) /(X[i0+1] - X[i0  ]);
      
      val = sm * Y[i0-1] + s0 * Y[i0] + s1 * Y[i0+1];
    }
  else if (order == 1)
    {
      double sm = (2.*x - X[i0  ] - X[i0+1]) /(X[i0-1] - X[i0  ]) /(X[i0-1] - X[i0+1]);
      double s0 = (2.*x - X[i0-1] - X[i0+1]) /(X[i0  ] - X[i0-1]) /(X[i0  ] - X[i0+1]);
      double s1 = (2.*x - X[i0-1] - X[i0  ]) /(X[i0+1] - X[i0-1]) /(X[i0+1] - X[i0  ]);
  
      val = sm * Y[i0-1] + s0 * Y[i0] + s1 * Y[i0+1];
    }
   else
    {
      printf ("FLUX::Interpolate: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ####################################################################
// 1D interpolation function with nonuniform grid and periodic function
// order = 0: Y(x)
// order = 1: dY/dx
// ####################################################################
double Flux::InterpolatePeriodic (int I, double* X, double* Y, double x, int order)
{  
  int i0 = 0, im, ip;
  for (int i = 1; i < I; i++)
    if (x > X[i])
      i0 = i;
  if (i0 < 0 || i0 > I-1)
    {
      printf ("FLUX::InterpolatePeriodic: Interpolation error: I = %3d i0 = %3d x = %11.4e X[0] = %11.4e X[I-1] = %11.4e\n",
	      I, i0, x, X[0], X[I-1]);
      exit (1);
    }
  if (x - X[i0] > 0.5 * (X[i0+1] - X[i0]))
    i0 += 1;
  im = i0 - 1;
  ip = i0 + 1;
  if (im == -1)
    im = I-1;
  if (ip == I)
    ip = 0;

  double val;
  if (order == 0)
    {
      double sm = (x - X[i0]) * (x - X[ip]) /(X[im] - X[i0]) /(X[im] - X[ip]);
      double s0 = (x - X[im]) * (x - X[ip]) /(X[i0] - X[im]) /(X[i0] - X[ip]);
      double s1 = (x - X[im]) * (x - X[i0]) /(X[ip] - X[im]) /(X[ip] - X[i0]);
      
      val = sm * Y[im] + s0 * Y[i0] + s1 * Y[ip];
    }
  else if (order == 1)
    {
      double sm = (2.*x - X[i0] - X[ip]) /(X[im] - X[i0]) /(X[im] - X[ip]);
      double s0 = (2.*x - X[im] - X[ip]) /(X[i0] - X[im]) /(X[i0] - X[ip]);
      double s1 = (2.*x - X[im] - X[i0]) /(X[ip] - X[im]) /(X[ip] - X[i0]);
  
      val = sm * Y[im] + s0 * Y[i0] + s1 * Y[ip];
    }
   else
    {
      printf ("FLUX::InterpolatePeriodic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ##############################################
// Function to interpolate Psi on uniform 2D grid
// order = 0: Psi
// order = 1: Psi_x
// order = 2: Psi_y
// order = 3: Psi_xx
// order = 4: Psi_yy
// ##############################################
double Flux::InterpolatePsi (double RR, double ZZ, int order)
{
  double i  = double (NRPTS-1) * (RR - RPTS[0]) /(RPTS[NRPTS-1] - RPTS[0]);
  int    i0 = int (i); 
  if (i0 < 0 || i0 > NRPTS-1)
    {
      printf ("FLUX::InterpolatePsi - Interpolation error: NRPTS = %3d i0 = %3d r = %11.4e RPTS[0] = %11.4e RPTS[NRPTS-1] = %11.4e\n", 
	      NRPTS, i0, RR, RPTS[0], RPTS[NRPTS-1]);
      exit (1);
    }
  if (i - double (i0) > 0.5)
    i0 += 1;
  if (i0 == 0)
    i0 += 1;
  if (i0 == NRPTS-1)
    i0 -= 1;

  double j  = double (NZPTS-1) * (ZZ - ZPTS[0]) /(ZPTS[NZPTS-1] - ZPTS[0]);
  int    j0 = int (j); 
  if (j0 < 0 || j0 > NZPTS-1)
    {
      printf ("FLUX::InterpolatePsi - Interpolation error: NZPTS = %3d j0 = %3d z = %11.4e ZPTS[0] = %11.4e ZPTS[NZPTS-1] = %11.4e\n", 
	      NZPTS, j0, ZZ, ZPTS[0], ZPTS[NZPTS-1]);
      exit (1);
    }
  if (j - double (j0) > 0.5)
    j0 += 1;
  if (j0 == 0)
    j0 += 1;
  if (j0 == NZPTS-1)
    j0 -= 1;

  double dR = (RPTS[1] - RPTS[0]);
  double dZ = (ZPTS[1] - ZPTS[0]);
  double x  = (RR - RPTS[i0]) /dR;
  double z  = (ZZ - ZPTS[j0]) /dZ;

  double val;
  if (order == 0)
    {
      double xm = + 0.5 * x * (x-1.);
      double x0 = - (x+1.) * (x-1.);
      double x1 = + 0.5 * x * (x+1.);
      double zm = + 0.5 * z * (z-1.);
      double z0 = - (z+1.) * (z-1.);
      double z1 = + 0.5 * z * (z+1.);

      val =
	+ xm * (+ zm * gsl_matrix_get (PSIARRAY, i0-1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0-1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (PSIARRAY, i0  , j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0  , j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (PSIARRAY, i0+1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0+1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0+1, j0+1));
    }
  else if (order == 1)
    {
      double xm = + x - 0.5;
      double x0 = - 2. * x;
      double x1 = + x + 0.5;
      double zm = + 0.5 * z * (z-1.);
      double z0 = - (z+1.) * (z-1.);
      double z1 = + 0.5 * z * (z+1.);
 
      val = 
	+ xm * (+ zm * gsl_matrix_get (PSIARRAY, i0-1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0-1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (PSIARRAY, i0  , j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0  , j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (PSIARRAY, i0+1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0+1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0+1, j0+1));

      val /= dR;
    }
  else if (order == 2)
    {
      double xm = + 0.5 * x * (x-1.);
      double x0 = - (x+1.) * (x-1.);
      double x1 = + 0.5 * x * (x+1.);
      double zm = + z - 0.5;
      double z0 = - 2. * z;
      double z1 = + z + 0.5;

      val =
	+ xm * (+ zm * gsl_matrix_get (PSIARRAY, i0-1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0-1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (PSIARRAY, i0  , j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0  , j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (PSIARRAY, i0+1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0+1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0+1, j0+1));

      val /= dZ;
    }
  else if (order == 3)
    {
      double xm = + 1.;
      double x0 = - 2.;
      double x1 = + 1.;
      double zm = + 0.5 * z * (z-1.);
      double z0 = - (z+1.) * (z-1.);
      double z1 = + 0.5 * z * (z+1.);
 
      val = 
	+ xm * (+ zm * gsl_matrix_get (PSIARRAY, i0-1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0-1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (PSIARRAY, i0  , j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0  , j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (PSIARRAY, i0+1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0+1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0+1, j0+1));
      
      val /= (dR*dR);
    }
  else if (order == 4)
    {
      double xm = + 0.5 * x * (x-1.);
      double x0 = - (x+1.) * (x-1.);
      double x1 = + 0.5 * x * (x+1.);
      double zm = + 1.;
      double z0 = - 2.;
      double z1 = + 1.;

      val =
	+ xm * (+ zm * gsl_matrix_get (PSIARRAY, i0-1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0-1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (PSIARRAY, i0  , j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0  , j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (PSIARRAY, i0+1, j0-1)
		+ z0 * gsl_matrix_get (PSIARRAY, i0+1, j0  )
		+ z1 * gsl_matrix_get (PSIARRAY, i0+1, j0+1));

      val /= (dZ*dZ);
    }
   else
    {
      printf ("FLUX::InterpolatePsi: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ###################################
// Function to evaluate dPsi/dR (R, Z)
// ###################################
double Flux::GetPsiR (double r, double z)
{
  return InterpolatePsi (r, z, 1);
}

// ###################################
// Function to evaluate dPsi/dZ (R, Z)
// ###################################
double Flux::GetPsiZ (double r, double z)
{
  return InterpolatePsi (r, z, 2);
}

// #######################################
// Function to evaluate d^2Psi/dR^2 (R, Z)
// #######################################
double Flux::GetPsiRR (double r, double z)
{
  return InterpolatePsi (r, z, 3);
}

// #######################################
// Function to evaluate d^2Psi/dZ^2 (R, Z)
// #######################################
double Flux::GetPsiZZ (double r, double z)
{
  return InterpolatePsi (r, z, 4);
}
