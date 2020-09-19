#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Interpolate        (int I, double* X, double* Y, double x, int order);
double InterpolateCubic   (double* X, double* Y, double x, int i0, int i1, int i2, int order);
double InterpolateQuartic (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

int main ()
{
  FILE* file = fopen ("../../Outputs/Stage3/profiles.txt", "r");

  int    I = 512;
  double c;

  double* psi = new double [I];
  double* r   = new double [I];
  double* ne  = new double [I];
  double* te  = new double [I];
  double* ti  = new double [I];
  double* ni  = new double [I];
  double* we  = new double [I];
  double* nn  = new double [I];

  double nna = 3.e17;
  double ln  = 1.2e-2;
  double a   = 0.855;

  for (int i = 0; i < I; i++)
    fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	    &psi[i], &r[i], &ne[i], &c, &te[i], &c, &c, &c, &ti[i], &c, &ni[i], &c, &c, &c, &we[i], &c,
	    &c, &c, &c, &c);
  
  fclose (file);

  printf ("\nProfile data:\n\n");
  
  for (int i = 0; i < I; i++)
    {
      nn [i] = nna * exp ((r[i]-1.)*a/ln) /1.e19;

      printf ("PsiN = %11.4e  n_e = %11.4e  T_e = %11.4e  T_i = %11.4e  n_I = %11.4e  w_E = %11.4e  n_n = %11.4e\n",
	      psi[i], ne[i], te[i], ti[i], ni[i], we[i], nn[i]);        
    }

  file = fopen ("145380_2500_uncertainty.txt", "r");

  int J = 401;

  double* psi1 = new double [J];
  double* ne1  = new double [J];
  double* te1  = new double [J];
  double* ti1  = new double [J];
  double* ni1  = new double [J];
  double* we1  = new double [J];

  for (int j = 0; j < J; j++)
    fscanf (file, "%lf %lf %lf %lf %lf %lf",
	    &psi1[j],  &ne1[j], &te1[j], &ni1[j], &ti1[j], &we1[j]);
  
  fclose (file);

  for (int j = 0; j < J; j++)
    {
      ne1[j] /= 1.e19;
      te1[j] /= 1.e3;
      ni1[j] /= 1.e19;
      ti1[j] /= 1.e3;
      we1[j] /= 1.e3;
    }
  
  printf ("\nError data:\n\n");
  
  for (int j = 0; j < J; j++)
    printf ("PsiN = %11.4e  n_e = %11.4e  T_e = %11.4e  T_i = %11.4e  n_I = %11.4e  w_E = %11.4e\n",
	    psi1[j], ne1[j], te1[j], ti1[j], ni1[j], we1[j]);

  double* ne_1 = new double [I];
  double* te_1 = new double [I];
  double* ti_1 = new double [I];
  double* ni_1 = new double [I];
  double* we_1 = new double [I];

  double* ne_2 = new double [I];
  double* te_2 = new double [I];
  double* ti_2 = new double [I];
  double* ni_2 = new double [I];
  double* we_2 = new double [I];

  for (int i = 0; i < I; i++)
    {
      ne_1[i] = ne[i] - 0.5 * Interpolate (J, psi1, ne1, psi[i], 0);
      ne_2[i] = ne[i] + 0.5 * Interpolate (J, psi1, ne1, psi[i], 0);
      te_1[i] = te[i] - 0.5 * Interpolate (J, psi1, te1, psi[i], 0);
      te_2[i] = te[i] + 0.5 * Interpolate (J, psi1, te1, psi[i], 0);
      ti_1[i] = ti[i] - 0.5 * Interpolate (J, psi1, ti1, psi[i], 0);
      ti_2[i] = ti[i] + 0.5 * Interpolate (J, psi1, ti1, psi[i], 0);
      ni_1[i] = ni[i] - 0.5 * Interpolate (J, psi1, ni1, psi[i], 0);
      ni_2[i] = ni[i] + 0.5 * Interpolate (J, psi1, ni1, psi[i], 0);
      we_1[i] = we[i] - 0.5 * Interpolate (J, psi1, we1, psi[i], 0);
      we_2[i] = we[i] + 0.5 * Interpolate (J, psi1, we1, psi[i], 0);
    }

  file = fopen ("profiles.out", "w");
  for (int i = 0; i < I; i++)
    fprintf (file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
	     psi[i], ne_1[i], ne[i], ne_2[i], te_1[i], te[i], te_2[i], ti_1[i], ti[i], ti_2[i],
	     ni_1[i], ni[i], ni_2[i], we_1[i], we[i], we_2[i], nn[i]);
  fclose (file);

  delete[] psi;  delete[] ne;  delete[] te;  delete[] ti;  delete[] ni;  delete[] we;
  delete[] psi1; delete[] ne1; delete[] te1; delete[] ti1; delete[] ni1; delete[] we1;

  delete[] ne_1; delete[] te_1; delete[] ti_1; delete[] ni_1; delete[] we_1;
  delete[] ne_2; delete[] te_2; delete[] ti_2; delete[] ni_2; delete[] we_2;
}

// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// order = 0: Y(x)
// order = 1: dY/dx
// ########################################################
double Interpolate (int I, double* X, double* Y, double x, int order)
{
  int index, cntrl;
  
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
  
  return val;
}

double InterpolateCubic (double* X, double* Y, double x, int i0, int i1, int i2, int order)
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
   else
    {
      printf ("FLUX::InterpolateCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double InterpolateQuartic (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
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
   else
    {
      printf ("FLUX::InterpolateQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}
