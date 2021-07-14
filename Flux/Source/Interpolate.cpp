// Interpolate.cpp

// PROGRAM ORGANIZATION:
//
// double Flux:: Interpolate                  (int I, double* X, double* Y, double x, int order)
// double Flux:: InterpolateCubic             (double* X, double* Y, double x, int i0, int i1, int i2, int order)
// double Flux:: InterpolateQuartic           (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
// double Flux:: InterpolatePeriodic          (int I, double* X, double* Y, double x, int order)
// double Flux:: InterpolatePeriodicQuartic   (int I, double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
// double Flux:: InterpolatePsi               (double RR, double ZZ, int order)
// double Flux:: InterpolatePsiCubicCubic     (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2, int order)
// double Flux:: InterpolatePsiQuarticCubic   (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int order)
// double Flux:: InterpolatePsiCubicQuartic   (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2, int j3, int order)
// double Flux:: InterpolatePsiQuarticQuartic (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order)
// double Flux:: GetPsiR                      (double r, double z)
// double Flux:: GetPsiZ                      (double r, double z)
// double Flux:: GetPsiRR                     (double r, double z)
// double Flux:: GetPsiRR                     (double r, double z)
// double Flux:: GetPsiRZ                     (double r, double z)
// double Flux:: GetPsiZZ                     (double r, double z)
// double Flux:: InterpolateQ                 (int I, double* X, double* Y, double x, int order, int method)
// double Flux:: InterpolateQQuartic          (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
// double Flux:: InterpolateQCubic            (double* X, double* Y, double x, int i0, int i1, int i2, int order)
// double Flux:: InterpolateQQuadratic        (double* X, double* Y, double x, int i0, int i1, int order)

#include "Flux.h"

// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// order = 0: Y(x)
// order = 1: dY/dx
// order = 2: d^2Y/dx^2
// ########################################################
double Flux::Interpolate (int I, double* X, double* Y, double x, int order)
{
  int index, cntrl = 0;
  
  if (x < X[0])
    {
      index = 0;
      cntrl = 2;
    }
  else if (x >= X[I-1])
    {
      index = I - 2;
      cntrl = 3;
    }
  else
    {
      for (int i = 0; i < I-1; i++)
	if (x >= X[i] && x < X[i+1])
	  {
	    index = i;

	    if (index == 0)
	      cntrl = 2;
	    else if (index == I-2)
	      cntrl = 3;
	    else
	      cntrl = 1;
	  }
    }

  double val;
  if (cntrl == 1)
    val = InterpolateQuartic (X, Y, x, index-1, index, index+1, index+2, order);
  else if (cntrl == 2)
    val = InterpolateCubic   (X, Y, x,          index, index+1, index+2, order);
  else if (cntrl == 3)
    val = InterpolateCubic   (X, Y, x, index-1, index, index+1,          order);
  else
    {
      printf ("FLUX::Interpolate: Error - cntrl = %1d\n", cntrl);
      exit (1);
    }
  
  return val;
}

double Flux::InterpolateCubic (double* X, double* Y, double x, int i0, int i1, int i2, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = (x - X[i0]) * (x - X[i2]) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = (x - X[i0]) * (x - X[i1]) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i1]) + (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = ((x - X[i0]) + (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = ((x - X[i0]) + (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else if (order == 2)
    {
      double s0 = 2. /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = 2. /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = 2. /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else
    {
      printf ("FLUX::InterpolateCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQuartic (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) * (x - X[i3]) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = (x - X[i0]) * (x - X[i2]) * (x - X[i3]) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = (x - X[i0]) * (x - X[i1]) * (x - X[i3]) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = (x - X[i0]) * (x - X[i1]) * (x - X[i2]) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i2]) * (x - X[i3]) + (x - X[i1]) * (x - X[i3]) + (x - X[i1]) * (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = ((x - X[i2]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = ((x - X[i1]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = ((x - X[i1]) * (x - X[i2]) + (x - X[i0]) * (x - X[i2]) + (x - X[i0]) * (x - X[i1])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 2)
    {
      double s0 = 2. * ((x - X[i1]) + (x - X[i2]) + (x - X[i3])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = 2. * ((x - X[i0]) + (x - X[i2]) + (x - X[i3])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i3])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i2])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else
    {
      printf ("FLUX::InterpolateQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ##############################################################################
// 1D interpolation function with nonuniform monotonic grid and periodic function
// Uses quartic interpolation
//
// order = 0: Y(x)
// order = 1: dY/dx
//
// ##############################################################################
double Flux::InterpolatePeriodic (int I, double* X, double* Y, double x, int order)
{
  int index;

  if (x < X[0])
      index = 0;
  else if (x >= X[I-1])
      index = I - 2;
  else
    {
      for (int i = 0; i < I-1; i++)
        if (x >= X[i] && x < X[i+1])
	    index = i;
    }
   
  double val = InterpolatePeriodicQuartic (I, X, Y, x, index-1, index, index+1, index+2, order);
 
  return val;
}

double Flux::InterpolatePeriodicQuartic (int I, double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    {
      if (i0 < 0)
	{
	  double Xi0 = X[0] - (X[1] - X[0]);

	  double s0 = (x - X[i1]) * (x - X[i2]) * (x - X[i3]) /(Xi0   - X[i1]) /(Xi0   - X[i2]) /(Xi0   - X[i3]);
	  double s1 = (x - Xi0  ) * (x - X[i2]) * (x - X[i3]) /(X[i1] - Xi0  ) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
	  double s2 = (x - Xi0  ) * (x - X[i1]) * (x - X[i3]) /(X[i2] - Xi0  ) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
	  double s3 = (x - Xi0  ) * (x - X[i1]) * (x - X[i2]) /(X[i3] - Xi0  ) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
	  
	  val = s0 * Y[I-1] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
	}
      else if (i3 > I-1)
	{
	  double Xi3 = X[I-1] + (X[I-1] - X[I-2]);

	  double s0 = (x - X[i1]) * (x - X[i2]) * (x - Xi3  ) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - Xi3  );
	  double s1 = (x - X[i0]) * (x - X[i2]) * (x - Xi3  ) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - Xi3  );
	  double s2 = (x - X[i0]) * (x - X[i1]) * (x - Xi3  ) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - Xi3  );
	  double s3 = (x - X[i0]) * (x - X[i1]) * (x - X[i2]) /(Xi3   - X[i0]) /(Xi3   - X[i1]) /(Xi3   - X[i2]);
	  
	  val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[0];
	}
      else
	{
	  double s0 = (x - X[i1]) * (x - X[i2]) * (x - X[i3]) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
	  double s1 = (x - X[i0]) * (x - X[i2]) * (x - X[i3]) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
	  double s2 = (x - X[i0]) * (x - X[i1]) * (x - X[i3]) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
	  double s3 = (x - X[i0]) * (x - X[i1]) * (x - X[i2]) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
	  
	  val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
	}
    }
  else if (order == 1)
    {
      if (i0 < 0)
	{
	  double Xi0 = X[0] - (X[1] - X[0]);
	  
	  double s0 = ((x - X[i2]) * (x - X[i3]) + (x - X[i1]) * (x - X[i3]) + (x - X[i1]) * (x - X[i2])) /(Xi0   - X[i1]) /(Xi0   - X[i2]) /(Xi0   - X[i3]);
	  double s1 = ((x - X[i2]) * (x - X[i3]) + (x - Xi0  ) * (x - X[i3]) + (x - Xi0  ) * (x - X[i2])) /(X[i1] - Xi0  ) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
	  double s2 = ((x - X[i1]) * (x - X[i3]) + (x - Xi0  ) * (x - X[i3]) + (x - Xi0  ) * (x - X[i1])) /(X[i2] - Xi0  ) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
	  double s3 = ((x - X[i1]) * (x - X[i2]) + (x - Xi0  ) * (x - X[i2]) + (x - Xi0  ) * (x - X[i1])) /(X[i3] - Xi0  ) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
	  
	  val = s0 * Y[I-1] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
	}
      else if (i3 > I-1)
	{
	  double Xi3 = X[I-1] + (X[I-1] - X[I-2]);
	   
	  double s0 = ((x - X[i2]) * (x - Xi3  ) + (x - X[i1]) * (x - Xi3  ) + (x - X[i1]) * (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - Xi3  );
	  double s1 = ((x - X[i2]) * (x - Xi3  ) + (x - X[i0]) * (x - Xi3  ) + (x - X[i0]) * (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - Xi3  );
	  double s2 = ((x - X[i1]) * (x - Xi3  ) + (x - X[i0]) * (x - Xi3  ) + (x - X[i0]) * (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - Xi3  );
	  double s3 = ((x - X[i1]) * (x - X[i2]) + (x - X[i0]) * (x - X[i2]) + (x - X[i0]) * (x - X[i1])) /(Xi3   - X[i0]) /(Xi3   - X[i1]) /(Xi3   - X[i2]);
	  
	  val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[0];
	}
      else
	{
	  double s0 = ((x - X[i2]) * (x - X[i3]) + (x - X[i1]) * (x - X[i3]) + (x - X[i1]) * (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
	  double s1 = ((x - X[i2]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
	  double s2 = ((x - X[i1]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
	  double s3 = ((x - X[i1]) * (x - X[i2]) + (x - X[i0]) * (x - X[i2]) + (x - X[i0]) * (x - X[i1])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
	  
	  val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
	}
    }
   else
    {
      printf ("FLUX::InterpolatePeriodicQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ###########################################################
// Function to interpolate Psi on uniform 2D grid
// Defaults to quartic interpolation
// Uses cubic/quardatic iterpolation close to  grid boundaries
//
// order = 0: Psi
// order = 1: Psi_r
// order = 2: Psi_z
// order = 3; Psi_rr
// order = 4; Psi_rz
// order = 5; Psi_zz
//
// ###########################################################
double Flux::InterpolatePsi (double RR, double ZZ, int order)
{
  int i0, ic = 0;
  
  if (RR < RPTS[0])
    {
      i0 = 0;
      ic = 2;
    }
  else if (RR >= RPTS[NRPTS-1])
    {
      i0 = NRPTS - 2;
      ic = 3;
    }
  else
    {
      for (int i = 0; i < NRPTS-1; i++)
	if (RR >= RPTS[i] && RR < RPTS[i+1])
	  {
	    i0 = i;

	    if (i0 == 0)
	      ic = 2;
	    else if (i0 == NRPTS-2)
	      ic = 3;
	    else
	      ic = 1;
	  }
    }

  int j0, jc = 0;
  
  if (ZZ < ZPTS[0])
    {
      j0 = 0;
      jc = 2;
    }
  else if (ZZ >= ZPTS[NZPTS-1])
    {
      j0 = NZPTS - 2;
      jc = 3;
    }
  else
    {
      for (int i = 0; i < NZPTS-1; i++)
	if (ZZ >= ZPTS[i] && ZZ < ZPTS[i+1])
	  {
	    j0 = i;

	    if (j0 == 0)
	      jc = 2;
	    else if (j0 == NZPTS-2)
	      jc = 3;
	    else
	      jc = 1;
	  }
    }

  double val;
  if      (ic == 1 && jc == 1)
    val = InterpolatePsiQuarticQuartic (RR, ZZ, i0-1, i0, i0+1, i0+2, j0-1, j0, j0+1, j0+2, order);
  else if (ic == 1 && jc == 2)
    val = InterpolatePsiQuarticCubic   (RR, ZZ, i0-1, i0, i0+1, i0+2,       j0, j0+1, j0+2, order);
  else if (ic == 1 && jc == 3)
    val = InterpolatePsiQuarticCubic   (RR, ZZ, i0-1, i0, i0+1, i0+2, j0-1, j0, j0+1,       order);
  else if (ic == 2 && jc == 1)
    val = InterpolatePsiCubicQuartic   (RR, ZZ,       i0, i0+1, i0+2, j0-1, j0, j0+1, j0+2, order);
  else if (ic == 3 && jc == 1)
    val = InterpolatePsiCubicQuartic   (RR, ZZ, i0-1, i0, i0+1,       j0-1, j0, j0+1, j0+2, order);
  else if (ic == 2 && jc == 2)
    val = InterpolatePsiCubicCubic     (RR, ZZ,       i0, i0+1, i0+2,       j0, j0+1, j0+2, order);
  else if (ic == 2 && jc == 3)
    val = InterpolatePsiCubicCubic     (RR, ZZ,       i0, i0+1, i0+2, j0-1, j0, j0+1,       order);
  else if (ic == 3 && jc == 2)
    val = InterpolatePsiCubicCubic     (RR, ZZ, i0-1, i0, i0+1,             j0, j0+1, j0+2, order);
  else if (ic == 3 && jc == 3)
    val = InterpolatePsiCubicCubic     (RR, ZZ, i0-1, i0, i0+1,       j0-1, j0, j0+1,       order);
  else
    {
      printf ("FLUX::InterpolatePsi: Error - ic = %1d, jc = %1d\n", ic, jc);
      exit (1);
    }	
  
  return val;
}

double Flux::InterpolatePsiCubicCubic (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2, int order)
{
  double val;

  double RR0 = RPTS[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;

  double ZZ0 = ZPTS[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);
      
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);
       
      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 3)
    {
      double r0 = 2. /(2.*dR2);
      double r1 = 2. /(-  dR2);
      double r2 = 2. /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);

      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);
       
      double z0 = 2. /(2.*dZ2);
      double z1 = 2. /(-  dZ2);
      double z2 = 2. /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else
    {
      printf ("FLUX::InterpolatePsiCubicCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolatePsiQuarticCubic (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int order)
{
  double val;
  
  double RR0 = RPTS[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;
  double RR3 = RR2 + dR;

  double ZZ0 = ZPTS[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
      
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);
      double val3 = z0 * PSIARRAY (i3, j0) + z1 * PSIARRAY (i3, j1) + z2 * PSIARRAY (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);
 
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);
      double val3 = z0 * PSIARRAY (i3, j0) + z1 * PSIARRAY (i3, j1) + z2 * PSIARRAY (i3, j2);
      
      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
       
      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);
      double val3 = z0 * PSIARRAY (i3, j0) + z1 * PSIARRAY (i3, j1) + z2 * PSIARRAY (i3, j2);
   
      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
   else if (order == 3)
    {
      double r0 = 2. * ((RR - RR1) + (RR - RR2) + (RR - RR3)) /(-6.*dR3);
      double r1 = 2. * ((RR - RR0) + (RR - RR2) + (RR - RR3)) /(+2.*dR3);
      double r2 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR3)) /(-2.*dR3);
      double r3 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR2)) /(+6.*dR3);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);
 
      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);
      double val3 = z0 * PSIARRAY (i3, j0) + z1 * PSIARRAY (i3, j1) + z2 * PSIARRAY (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);

      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);
      double val3 = z0 * PSIARRAY (i3, j0) + z1 * PSIARRAY (i3, j1) + z2 * PSIARRAY (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val + r3 * val3;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
       
      double z0 = 2. /(2.*dZ2);
      double z1 = 2. /(-  dZ2);
      double z2 = 2. /(2.*dZ2);

      double val0 = z0 * PSIARRAY (i0, j0) + z1 * PSIARRAY (i0, j1) + z2 * PSIARRAY (i0, j2);
      double val1 = z0 * PSIARRAY (i1, j0) + z1 * PSIARRAY (i1, j1) + z2 * PSIARRAY (i1, j2);
      double val2 = z0 * PSIARRAY (i2, j0) + z1 * PSIARRAY (i2, j1) + z2 * PSIARRAY (i2, j2);
      double val3 = z0 * PSIARRAY (i3, j0) + z1 * PSIARRAY (i3, j1) + z2 * PSIARRAY (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else
    {
      printf ("FLUX::InterpolatePsiQuarticCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolatePsiCubicQuartic (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2, int j3, int order)
{
  double val;

  double RR0 = RPTS[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;

  double ZZ0 = ZPTS[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;
  double ZZ3 = ZZ2 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);
 
      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);

      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 3)
    {
      double r0 = 2. /(2.*dR2);
      double r1 = 2. /(-  dR2);
      double r2 = 2. /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);
      
      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);

      double z0 = 2. * ((ZZ - ZZ1) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(-6.*dZ3);
      double z1 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(+2.*dZ3);
      double z2 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ3)) /(-2.*dZ3);
      double z3 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ2)) /(+6.*dZ3);
 
      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else
    {
      printf ("FLUX::InterpolatePsiCubicQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolatePsiQuarticQuartic (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order)
{
  double val;

  double RR0 = RPTS[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;
  double RR3 = RR2 + dR;

  double ZZ0 = ZPTS[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;
  double ZZ3 = ZZ2 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0) + r3 * PSIARRAY (i3, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1) + r3 * PSIARRAY (i3, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2) + r3 * PSIARRAY (i3, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3) + r3 * PSIARRAY (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0) + r3 * PSIARRAY (i3, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1) + r3 * PSIARRAY (i3, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2) + r3 * PSIARRAY (i3, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3) + r3 * PSIARRAY (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
      
      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0) + r3 * PSIARRAY (i3, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1) + r3 * PSIARRAY (i3, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2) + r3 * PSIARRAY (i3, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3) + r3 * PSIARRAY (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 3)
    {
      double r0 = 2. * ((RR - RR1) + (RR - RR2) + (RR - RR3)) /(-6.*dR3);
      double r1 = 2. * ((RR - RR0) + (RR - RR2) + (RR - RR3)) /(+2.*dR3);
      double r2 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR3)) /(-2.*dR3);
      double r3 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR2)) /(+6.*dR3);
  
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0) + r3 * PSIARRAY (i3, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1) + r3 * PSIARRAY (i3, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2) + r3 * PSIARRAY (i3, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3) + r3 * PSIARRAY (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);

      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0) + r3 * PSIARRAY (i3, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1) + r3 * PSIARRAY (i3, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2) + r3 * PSIARRAY (i3, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3) + r3 * PSIARRAY (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);

      double z0 = 2. * ((ZZ - ZZ1) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(-6.*dZ3);
      double z1 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(+2.*dZ3);
      double z2 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ3)) /(-2.*dZ3);
      double z3 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ2)) /(+6.*dZ3);
 
      double val0 = r0 * PSIARRAY (i0, j0) + r1 * PSIARRAY (i1, j0) + r2 * PSIARRAY (i2, j0) + r3 * PSIARRAY (i3, j0);
      double val1 = r0 * PSIARRAY (i0, j1) + r1 * PSIARRAY (i1, j1) + r2 * PSIARRAY (i2, j1) + r3 * PSIARRAY (i3, j1);
      double val2 = r0 * PSIARRAY (i0, j2) + r1 * PSIARRAY (i1, j2) + r2 * PSIARRAY (i2, j2) + r3 * PSIARRAY (i3, j2);
      double val3 = r0 * PSIARRAY (i0, j3) + r1 * PSIARRAY (i1, j3) + r2 * PSIARRAY (i2, j3) + r3 * PSIARRAY (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else
    {
      printf ("FLUX::InterpolatePsiQuarticQuartic: Error - order = %1d\n", order);
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
// Function to evaluate d^2Psi/dRdZ (R, Z)
// #######################################
double Flux::GetPsiRZ (double r, double z)
{
  return InterpolatePsi (r, z, 4);
}

// #######################################
// Function to evaluate d^2Psi/dZ^2 (R, Z)
// #######################################
double Flux::GetPsiZZ (double r, double z)
{
  return InterpolatePsi (r, z, 5);
}

// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// specifically written to interpolate and extrapolate Q
// profile
//
// order = 0: Y(x)
// order = 1: dY/dx
// order = 2: d^2Y/dx^2
// order = 3: d^3Y/dx^3
//
// method = 2: quadratic (linear) interpolation
// method = 3: cubic interpolation
// method = 4: quartic interpolation
//
// ########################################################
double Flux::InterpolateQ (int I, double* X, double* Y, double x, int order, int method)
{
  int    index, cntrl = 0;
  double val;

  if (method == 4)
    {
      if (x < X[1])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (x >= X[I-2])
	{
	  index = I - 1;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < I-1; i++)
	    if (x >= X[i] && x < X[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index > I-3)
		  {
		    cntrl = 3;
		    index = I - 1;
		  }
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	val = InterpolateQQuartic (X, Y, x, index-1, index, index+1, index+2, order);
      else if (cntrl == 2)
	val = InterpolateQQuartic (X, Y, x, index, index+1, index+2, index+4, order);
      else if (cntrl == 3)
	val = InterpolateQQuartic (X, Y, x, index-3, index-2, index-1, index, order);
      else
	{
	  printf ("FLUX::InterpolateQ: Error - method = %1d, ctrl = %1d\n", method, cntrl);
	  exit (1);
	}	
    }
  else if (method == 3)
    {
       if (x < X[1])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (x >= X[I-1])
	{
	  index = I - 1;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < I-1; i++)
	    if (x >= X[i] && x < X[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index > I-2)
		  {
		    cntrl = 3;
		    index = I - 1;
		  }
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	val = InterpolateQCubic (X, Y, x, index-1, index, index+1, order);
      else if (cntrl == 2)
	val = InterpolateQCubic (X, Y, x, index, index+1, index+2, order);
      else if (cntrl == 3)
	val = InterpolateQCubic (X, Y, x, index-2, index-1, index, order);
      else
	{
	  printf ("FLUX::InterpolateQ: Error - method = %1d, ctrl = %1d\n", method, cntrl);
	  exit (1);
	}	
    }
  else if (method == 2)
    {
      if (x < X[0])
	 {
	   index = 0;
	   cntrl = 2;
	 }
      else if (x >= X[I-1])
	{
	  index = I - 1;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < I-1; i++)
	    if (x >= X[i] && x < X[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index > I-2)
		  {
		    cntrl = 3;
		    index = I - 1;
		  }
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	val = InterpolateQQuadratic (X, Y, x, index, index+1, order);
      else if (cntrl == 2)
	val = InterpolateQQuadratic (X, Y, x, index, index+1, order);
      else if (cntrl == 3)
	val = InterpolateQQuadratic (X, Y, x, index-1, index, order);
      else
	{
	  printf ("FLUX::InterpolateQ: Error - method = %1d, ctrl = %1d\n", method, cntrl);
	  exit (1);
	}
    }
  
  return val;
}

double Flux::InterpolateQQuartic (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) * (x - X[i3]) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = (x - X[i0]) * (x - X[i2]) * (x - X[i3]) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = (x - X[i0]) * (x - X[i1]) * (x - X[i3]) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = (x - X[i0]) * (x - X[i1]) * (x - X[i2]) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i2]) * (x - X[i3]) + (x - X[i1]) * (x - X[i3]) + (x - X[i1]) * (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = ((x - X[i2]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = ((x - X[i1]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = ((x - X[i1]) * (x - X[i2]) + (x - X[i0]) * (x - X[i2]) + (x - X[i0]) * (x - X[i1])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 2)
    {
      double s0 = 2. * ((x - X[i1]) + (x - X[i2]) + (x - X[i3])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = 2. * ((x - X[i0]) + (x - X[i2]) + (x - X[i3])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i3])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i2])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 3)
    {
      double s0 = 6. /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = 6. /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = 6. /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = 6. /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else
    {
      printf ("FLUX::InterpolateQQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQCubic (double* X, double* Y, double x, int i0, int i1, int i2, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = (x - X[i0]) * (x - X[i2]) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = (x - X[i0]) * (x - X[i1]) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i1]) + (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = ((x - X[i0]) + (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = ((x - X[i0]) + (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else if (order == 2)
    {
      double s0 = 2. /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = 2. /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = 2. /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else
    {
      printf ("FLUX::InterpolateQCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQQuadratic (double* X, double* Y, double x, int i0, int i1, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) /(X[i0] - X[i1]);
      double s1 = (x - X[i0]) /(X[i1] - X[i0]);
       
      val = s0 * Y[i0] + s1 * Y[i1];
    }
  else if (order == 1)
    {
      double s0 = 1. /(X[i0] - X[i1]);
      double s1 = 1. /(X[i1] - X[i0]);
  
      val = s0 * Y[i0] + s1 * Y[i1];
    }
  else
    {
      printf ("FLUX::InterpolateQQuadratic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}
