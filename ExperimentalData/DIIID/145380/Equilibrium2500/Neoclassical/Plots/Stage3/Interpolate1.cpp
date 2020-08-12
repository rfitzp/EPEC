#include <stdio.h>
#include <stdlib.h>

double Interpolate        (int I, double* X, double* Y, double x, int order);
double InterpolateCubic   (double* X, double* Y, double x, int i0, int i1, int i2, int order);
double InterpolateQuartic (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

int main ()
{
  FILE* file = fopen ("coeff.txt", "r");

  int I = 50;

  double* psi   = new double [I];
  double* chin  = new double [I];
  double* dchin = new double [I];
  double* chie  = new double [I];
  double* dchie = new double [I];
  double* chip  = new double [I];
  double* dchip = new double [I];

  for (int i = 0; i < I; i++)
    fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
	    &psi[i], &chin[i], &dchin[i], &chie[i], &dchie[i], &chip[i], &dchip[i]);
  
  fclose (file);

  printf ("\nProfile data:\n\n");
  
  for (int i = 0; i < I; i++)
    printf ("PsiN = %11.4e  chi_n = %11.4e  chi_e = %11.4e  chi_p = %11.4e\n",
	    psi[i], chin[i], chie[i], chip[i]);

  int J = 512;

  double* psi1    = new double [J];
  double* chin_1  = new double [J];
  double* dchin_1 = new double [J];
  double* chie_1  = new double [J];
  double* dchie_1 = new double [J];
  double* chip_1  = new double [J];
  double* dchip_1 = new double [J];

  for (int j = 0; j < J; j++)
    psi1[j] = double (j) /double (J-1);
  
  for (int j = 0; j < J; j++)
    {
      chin_1 [j] = Interpolate (I, psi, chin,  psi1[j], 0);
      dchin_1[j] = Interpolate (I, psi, dchin, psi1[j], 0);
      chie_1 [j] = Interpolate (I, psi, chie,  psi1[j], 0);
      dchie_1[j] = Interpolate (I, psi, dchie, psi1[j], 0);
      chip_1 [j] = Interpolate (I, psi, chip,  psi1[j], 0);
      dchip_1[j] = Interpolate (I, psi, dchip, psi1[j], 0);
    }

  file = fopen ("diffusivity.out", "w");
  for (int j = 0; j < J; j++)
    fprintf (file, "%e %e %e %e %e %e %e\n",
	     psi1[j], chin_1[j], dchin_1[j], chie_1[j], dchie_1[j], chip_1[j], dchip_1[j]);
  fclose (file);

  delete[] psi;  delete[] chin;   delete[] dchin;   delete[] chie;   delete[] dchie;   delete[] chip;   delete[] dchip;
  delete[] psi1; delete[] chin_1; delete[] dchin_1; delete[] chie_1; delete[] dchie_1; delete[] chip_1; delete[] dchip_1;
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
