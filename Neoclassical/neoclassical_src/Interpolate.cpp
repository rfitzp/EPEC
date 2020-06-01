// Interpolate.cpp

#include "Neoclassical.h"
#include "Field.h"

// ##############################################
// 1D interpolation function with nonuniform grid
// order = 0: Y(x)
// order = 1: dY/dx
// ##############################################
double Neoclassical::Interpolate (int I, Array<double,1> X, Array<double,1> Y, double x, int order)
{  
  int i0 = 0;
  for (int i = 1; i < I; i++)
    if (x > X (i))
      i0 = i;
  if (i0 < 0 || i0 > I-1)
    {
      printf ("NEOCLASSICAL::Interpolate: Interpolation error: I = %3d i0 = %3d x = %11.4e X(0) = %11.4e X(I-1) = %11.4e\n",
	      I, i0, x, X (0), X (I-1));
      exit (1);
    }
  if (x - X (i0) > 0.5 * (X (i0+1) - X (i0)))
    i0 += 1;
  if (i0 == 0)
    i0 += 1;
  if (i0 == I-1)
    i0 -= 1;

  double val;
  if (order == 0)
    {
      double sm = (x - X (i0  )) * (x - X (i0+1)) /(X (i0-1) - X (i0  )) /(X (i0-1) - X (i0+1));
      double s0 = (x - X (i0-1)) * (x - X (i0+1)) /(X (i0  ) - X (i0-1)) /(X (i0  ) - X (i0+1));
      double s1 = (x - X (i0-1)) * (x - X (i0  )) /(X (i0+1) - X (i0-1)) /(X (i0+1) - X (i0  ));
      
      val = sm * Y (i0-1) + s0 * Y (i0) + s1 * Y (i0+1);
    }
  else if (order == 1)
    {
      double sm = (2.*x - X (i0  ) - X (i0+1)) /(X (i0-1) - X (i0  )) /(X (i0-1) - X (i0+1));
      double s0 = (2.*x - X (i0-1) - X (i0+1)) /(X (i0  ) - X (i0-1)) /(X (i0  ) - X (i0+1));
      double s1 = (2.*x - X (i0-1) - X (i0  )) /(X (i0+1) - X (i0-1)) /(X (i0+1) - X (i0  ));
  
      val = sm * Y (i0-1) + s0 * Y (i0) + s1 * Y (i0+1);
    }
  else
    {
      printf ("NEOCLASSICAL::Interpolate: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ##############################################
// 1D interpolation function with nonuniform grid
// order = 0: Y(x)
// order = 1: dY/dx
// ##############################################
double Neoclassical::InterpolateField (Field& F, double x, int order)
{
  double I = F.GetN ();
  
  int i0 = 0;
  for (int i = 1; i < I; i++)
    if (x > F.GetX (i))
      i0 = i;
  if (i0 < 0 || i0 > I-1)
    {
      printf ("NEOCLASSICAL::InterpolateField - Interpolation error: I = %3d i0 = %3d x = %11.4e X(0) = %11.4e X(I-1) = %11.4e\n",
	      I, i0, x, F.GetX (0), F.GetX (I-1));
      exit (1);
    }
  if (x - F.GetX (i0) > 0.5 * (F.GetX (i0+1) - F.GetX (i0)))
    i0 += 1;
  if (i0 == 0)
    i0 += 1;
  if (i0 == I-1)
    i0 -= 1;

  double val;
  if (order == 0)
    {
      double sm = (x - F.GetX (i0  )) * (x - F.GetX (i0+1)) /(F.GetX (i0-1) - F.GetX (i0  )) /(F.GetX (i0-1) - F.GetX (i0+1));
      double s0 = (x - F.GetX (i0-1)) * (x - F.GetX (i0+1)) /(F.GetX (i0  ) - F.GetX (i0-1)) /(F.GetX (i0  ) - F.GetX (i0+1));
      double s1 = (x - F.GetX (i0-1)) * (x - F.GetX (i0  )) /(F.GetX (i0+1) - F.GetX (i0-1)) /(F.GetX (i0+1) - F.GetX (i0  ));
      
      val = sm * F.GetY (i0-1) + s0 * F.GetY (i0) + s1 * F.GetY (i0+1);
    }
  else if (order == 1)
    {
      double sm = (x - F.GetX (i0  )) * (x - F.GetX (i0+1)) /(F.GetX (i0-1) - F.GetX (i0  )) /(F.GetX (i0-1) - F.GetX (i0+1));
      double s0 = (x - F.GetX (i0-1)) * (x - F.GetX (i0+1)) /(F.GetX (i0  ) - F.GetX (i0-1)) /(F.GetX (i0  ) - F.GetX (i0+1));
      double s1 = (x - F.GetX (i0-1)) * (x - F.GetX (i0  )) /(F.GetX (i0+1) - F.GetX (i0-1)) /(F.GetX (i0+1) - F.GetX (i0  ));
  
      val = sm * F.GetdYdX (i0-1) + s0 * F.GetdYdX (i0) + s1 * F.GetdYdX (i0+1);
    }
   else
    {
      printf ("NEOCLASSICAL::InterpolateField: Error - order = %1d\n", order);
      exit (1);
    }


  return val;
}

