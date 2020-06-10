// Interpolate.cpp

// PROGRAM ORGANIZATION:

// double Neoclassical:: Interpolate             (int I, Array<double,1> X, Array<double,1> Y, double x, int order)
// double Neoclassical:: InterpolateCubic        (Array<double,1> X, Array<double,1> Y, double x, int i0, int i1, int i2, int order)
// double Neoclassical:: InterpolateQuartic      (Array<double,1> X, Array<double,1> Y, double x, int i0, int i1, int i2, int i3, int order)
// double Neoclassical:: InterpolateField        (Field& F, double x, int order)
// double Neoclassical:: InterpolateFieldCubic   (Field& F, double x, int i0, int i1, int i2, int order)
// double Neoclassical:: InterpolateFieldQuartic (Field& F, double x, int i0, int i1, int i2, int i3, int order)

#include "Neoclassical.h"

// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// order = 0: Y(x)
// order = 1: dY/dx
// ########################################################
double Neoclassical::Interpolate (int I, Array<double,1> X, Array<double,1> Y, double x, int order)
{
  int index, cntrl;
  
  if (x < X(0))
    {
      index = 0;
      cntrl = 2;
    }
  else if (x >= X(I-1))
    {
      index = I - 2;
      cntrl = 3;
    }
  else
    {
      for (int i = 0; i < I-1; i++)
	if (x >= X(i) && x < X(i+1))
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
  
  return val;
}

double Neoclassical::InterpolateCubic (Array<double,1> X, Array<double,1> Y, double x, int i0, int i1, int i2, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X(i1)) * (x - X(i2)) /(X(i0) - X(i1)) /(X(i0) - X(i2));
      double s1 = (x - X(i0)) * (x - X(i2)) /(X(i1) - X(i0)) /(X(i1) - X(i2));
      double s2 = (x - X(i0)) * (x - X(i1)) /(X(i2) - X(i0)) /(X(i2) - X(i1));
      
      val = s0 * Y(i0) + s1 * Y(i1) + s2 * Y(i2);
    }
  else if (order == 1)
    {
      double s0 = ((x - X(i1)) + (x - X(i2))) /(X(i0) - X(i1)) /(X(i0) - X(i2));
      double s1 = ((x - X(i0)) + (x - X(i2))) /(X(i1) - X(i0)) /(X(i1) - X(i2));
      double s2 = ((x - X(i0)) + (x - X(i1))) /(X(i2) - X(i0)) /(X(i2) - X(i1));
  
      val = s0 * Y(i0) + s1 * Y(i1) + s2 * Y(i2);
    }
   else
    {
      printf ("NEOCLASSICAL::InterpolateCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Neoclassical::InterpolateQuartic (Array<double,1> X, Array<double,1> Y, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X(i1)) * (x - X(i2)) * (x - X(i3)) /(X(i0) - X(i1)) /(X(i0) - X(i2)) /(X(i0) - X(i3));
      double s1 = (x - X(i0)) * (x - X(i2)) * (x - X(i3)) /(X(i1) - X(i0)) /(X(i1) - X(i2)) /(X(i1) - X(i3));
      double s2 = (x - X(i0)) * (x - X(i1)) * (x - X(i3)) /(X(i2) - X(i0)) /(X(i2) - X(i1)) /(X(i2) - X(i3));
      double s3 = (x - X(i0)) * (x - X(i1)) * (x - X(i2)) /(X(i3) - X(i0)) /(X(i3) - X(i1)) /(X(i3) - X(i2));
      
      val = s0 * Y(i0)+ s1 * Y(i1) + s2 * Y(i2) + s3 * Y(i3);
    }
  else if (order == 1)
    {
      double s0 = ((x - X(i2)) * (x - X(i3)) + (x - X(i1)) * (x - X(i3)) + (x - X(i1)) * (x - X(i2))) /(X(i0) - X(i1)) /(X(i0) - X(i2)) /(X(i0) - X(i3));
      double s1 = ((x - X(i2)) * (x - X(i3)) + (x - X(i0)) * (x - X(i3)) + (x - X(i0)) * (x - X(i2))) /(X(i1) - X(i0)) /(X(i1) - X(i2)) /(X(i1) - X(i3));
      double s2 = ((x - X(i1)) * (x - X(i3)) + (x - X(i0)) * (x - X(i3)) + (x - X(i0)) * (x - X(i1))) /(X(i2) - X(i0)) /(X(i2) - X(i1)) /(X(i2) - X(i3));
      double s3 = ((x - X(i1)) * (x - X(i2)) + (x - X(i0)) * (x - X(i2)) + (x - X(i0)) * (x - X(i1))) /(X(i3) - X(i0)) /(X(i3) - X(i1)) /(X(i3) - X(i2));
      
      val = s0 * Y(i0) + s1 * Y(i1) + s2 * Y(i2) + s3 * Y(i3);
    }
   else
    {
      printf ("NEOCLASSICAL::InterpolateQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}
    
// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// order = 0: Y(x)
// order = 1: dY/dx
// #######################################################
double Neoclassical::InterpolateField (Field& F, double x, int order)
{
  double I = F.GetN ();
  
  int index, cntrl;
  
  if (x < F.GetX(0))
    {
      index = 0;
      cntrl = 2;
    }
  else if (x >= F.GetX(I-1))
    {
      index = I - 2;
      cntrl = 3;
    }
  else
    {
      for (int i = 0; i < I-1; i++)
	if (x >= F.GetX(i) && x < F.GetX(i+1))
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
    val = InterpolateFieldQuartic (F, x, index-1, index, index+1, index+2, order);
  else if (cntrl == 2)
    val = InterpolateFieldCubic   (F, x,          index, index+1, index+2, order);
  else if (cntrl == 3)
    val = InterpolateFieldCubic   (F, x, index-1, index, index+1,          order);
  
  return val;
}

double Neoclassical::InterpolateFieldCubic (Field& F, double x, int i0, int i1, int i2, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - F.GetX(i1)) * (x - F.GetX(i2)) /(F.GetX(i0) - F.GetX(i1)) /(F.GetX(i0) - F.GetX(i2));
      double s1 = (x - F.GetX(i0)) * (x - F.GetX(i2)) /(F.GetX(i1) - F.GetX(i0)) /(F.GetX(i1) - F.GetX(i2));
      double s2 = (x - F.GetX(i0)) * (x - F.GetX(i1)) /(F.GetX(i2) - F.GetX(i0)) /(F.GetX(i2) - F.GetX(i1));
      
      val = s0 * F.GetY(i0) + s1 * F.GetY(i1) + s2 * F.GetY(i2);
    }
  else if (order == 1)
    {
      double s0 = (x - F.GetX(i1)) * (x - F.GetX(i2)) /(F.GetX(i0) - F.GetX(i1)) /(F.GetX(i0) - F.GetX(i2));
      double s1 = (x - F.GetX(i0)) * (x - F.GetX(i2)) /(F.GetX(i1) - F.GetX(i0)) /(F.GetX(i1) - F.GetX(i2));
      double s2 = (x - F.GetX(i0)) * (x - F.GetX(i1)) /(F.GetX(i2) - F.GetX(i0)) /(F.GetX(i2) - F.GetX(i1));
      
      val = s0 * F.GetdYdX(i0) + s1 * F.GetdYdX(i1) + s2 * F.GetdYdX(i2);
    }
   else
    {
      printf ("NEOCLASSICAL::InterpolateFieldCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Neoclassical::InterpolateFieldQuartic (Field& F, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    { 
      double s0 = (x - F.GetX(i1)) * (x - F.GetX(i2)) * (x - F.GetX(i3)) /(F.GetX(i0) - F.GetX(i1)) /(F.GetX(i0) - F.GetX(i2)) /(F.GetX(i0) - F.GetX(i3));
      double s1 = (x - F.GetX(i0)) * (x - F.GetX(i2)) * (x - F.GetX(i3)) /(F.GetX(i1) - F.GetX(i0)) /(F.GetX(i1) - F.GetX(i2)) /(F.GetX(i1) - F.GetX(i3));
      double s2 = (x - F.GetX(i0)) * (x - F.GetX(i1)) * (x - F.GetX(i3)) /(F.GetX(i2) - F.GetX(i0)) /(F.GetX(i2) - F.GetX(i1)) /(F.GetX(i2) - F.GetX(i3));
      double s3 = (x - F.GetX(i0)) * (x - F.GetX(i1)) * (x - F.GetX(i2)) /(F.GetX(i3) - F.GetX(i0)) /(F.GetX(i3) - F.GetX(i1)) /(F.GetX(i3) - F.GetX(i2));
      
      val = s0 * F.GetY(i0)+ s1 * F.GetY(i1) + s2 * F.GetY(i2) + s3 * F.GetY(i3);
    }
  else if (order == 1)
    {
      double s0 = (x - F.GetX(i1)) * (x - F.GetX(i2)) * (x - F.GetX(i3)) /(F.GetX(i0) - F.GetX(i1)) /(F.GetX(i0) - F.GetX(i2)) /(F.GetX(i0) - F.GetX(i3));
      double s1 = (x - F.GetX(i0)) * (x - F.GetX(i2)) * (x - F.GetX(i3)) /(F.GetX(i1) - F.GetX(i0)) /(F.GetX(i1) - F.GetX(i2)) /(F.GetX(i1) - F.GetX(i3));
      double s2 = (x - F.GetX(i0)) * (x - F.GetX(i1)) * (x - F.GetX(i3)) /(F.GetX(i2) - F.GetX(i0)) /(F.GetX(i2) - F.GetX(i1)) /(F.GetX(i2) - F.GetX(i3));
      double s3 = (x - F.GetX(i0)) * (x - F.GetX(i1)) * (x - F.GetX(i2)) /(F.GetX(i3) - F.GetX(i0)) /(F.GetX(i3) - F.GetX(i1)) /(F.GetX(i3) - F.GetX(i2));
      
      val = s0 * F.GetdYdX(i0)+ s1 * F.GetdYdX(i1) + s2 * F.GetdYdX(i2) + s3 * F.GetdYdX(i3);
    }
   else
    {
      printf ("NEOCLASSICAL::InterpolateFieldQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

