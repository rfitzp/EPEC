// pFileInterpolate.h

#include "Neoclassical.h"

// ##############################
// Function to interpolate Fields
// ##############################
void Neoclassical::FieldInterpolate (Field& Field1, Field& Field2, Field& Field3, Field& Field, double weight1, double weight2, double weight3)
{
  int n;
  int n1 = Field1.GetN ();
  int n2 = Field2.GetN ();
  int n3 = Field3.GetN ();

  if (n1 == n2 && n2 == n3)
    {
      n = n1;
    }
  else
    {
      printf ("Error - Field size mismatch\n");
      exit (1);
    }

  Field.resize (n);

  double x1, y1, dydx1;
  double x2, y2, dydx2;
  double x3, y3, dydx3;
  double x,  y,  dydx;
  for (int i = 0; i < n; i++)
    {
      Field1.PullData (i, x1, y1, dydx1);
      Field2.PullData (i, x2, y2, dydx2);
      Field3.PullData (i, x3, y3, dydx3);

      x    = weight1 * x1    + weight2 * x2    + weight3 * x3;
      y    = weight1 * y1    + weight2 * y2    + weight3 * y3;
      dydx = weight1 * dydx1 + weight2 * dydx2 + weight3 * dydx3;

      Field.PushData (i, x, y, dydx);
    }
}

// ###############################
// Functions to interpolate pFiles
// ###############################
void Neoclassical::pFileInterpolate (char* pFile1, double time1, char* pFile2, double time2, char* pFile3, double time3, char* pFile, double time)
{
  int    n;
  double x, y, dydx;
  char   s1[30], s2[30], s3[30], s4[30], s5[30], s6[30];

  // ################
  // Read first pFile
  // ################
 
  // Check for existence of pfile
  FILE* file = OpenFiler (pFile1);
 
  // Read ne field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ne_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ne_1\n");
	  exit (1);
	}
      else
	{
	  ne_1.PushData (i, x, y, dydx);
	}
    }

  // Read Te field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field Te_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading Te_1\n");
	  exit (1);
	}
      else
	{
	  Te_1.PushData (i, x, y, dydx);
	}
    }

  // Read ni field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ni_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ni_1\n");
	  exit (1);
	}
      else
	{
	  ni_1.PushData (i, x, y, dydx);
	}
    }

  // Read Ti field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field Ti_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading Ti_1\n");
	  exit (1);
	}
      else
	{
	  Ti_1.PushData (i, x, y, dydx);
	}
    }

  // Read nb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field nb_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading nb_1\n");
	  exit (1);
	}
      else
	{
	  nb_1.PushData (i, x, y, dydx);
	}
    }

  // Read pb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field pb_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading pb_1\n");
	  exit (1);
	}
      else
	{
	  pb_1.PushData (i, x, y, dydx);
	}
    }

  // Read ptot field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ptot_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ptot_1\n");
	  exit (1);
	}
      else
	{
	  ptot_1.PushData (i, x, y, dydx);
	}
    }

  // Read omeg field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omeg_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omeg_1\n");
	  exit (1);
	}
      else
	{
	  omeg_1.PushData (i, x, y, dydx);
	}
    }

  // Read omegp field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omegp_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omegp_1\n");
	  exit (1);
	}
      else
	{
	  omegp_1.PushData (i, x, y, dydx);
	}
    }

  // Read omgvb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omgvb_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgvb_1\n");
	  exit (1);
	}
      else
	{
	  omgvb_1.PushData (i, x, y, dydx);
	}
    }

  // Read omgpp field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omgpp_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgpp_1\n");
	  exit (1);
	}
      else
	{
	  omgpp_1.PushData (i, x, y, dydx);
	}
    }
  
  // Read omgeb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field wE_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgeb_1\n");
	  exit (1);
	}
      else
	{
	  wE_1.PushData (i, x, y, dydx);
	}
    }

  // Read er field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field er_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading er_1\n");
	  exit (1);
	}
      else
	{
	  er_1.PushData (i, x, y, dydx);
	}
    }

  // Read ommvb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ommvb_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ommvb_1\n");
	  exit (1);
	}
      else
	{
	  ommvb_1.PushData (i, x, y, dydx);
	}
    }

  // Read ommpp field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ommpp_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ommpp_1\n");
	  exit (1);
	}
      else
	{
	  ommpp_1.PushData (i, x, y, dydx);
	}
    }

  // Read omevb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omevb_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omevb_1\n");
	  exit (1);
	}
      else
	{
	  omevb_1.PushData (i, x, y, dydx);
	}
    }

  // Read omepp field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omepp_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omepp_1\n");
	  exit (1);
	}
      else
	{
	  omepp_1.PushData (i, x, y, dydx);
	}
    }

  // Read kpol field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field kpol_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading kpol_1\n");
	  exit (1);
	}
      else
	{
	  kpol_1.PushData (i, x, y, dydx);
	}
    }
  
  // Read omghb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omghb_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omghb_1\n");
	  exit (1);
	}
      else
	{
	  omghb_1.PushData (i, x, y, dydx);
	}
    }

  // Read nI field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field nI_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading nz1_1\n");
	  exit (1);
	}
      else
	{
	  nI_1.PushData (i, x, y, dydx);
	}
    }

  // Read vtor1 field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field vtor1_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading vtor1_1\n");
	  exit (1);
	}
      else
	{
	  vtor1_1.PushData (i, x, y, dydx);
	}
    }

  // Read vpol1 field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field vpol1_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading vpol1_1\n");
	  exit (1);
	}
      else
	{
	  vpol1_1.PushData (i, x, y, dydx);
	}
    }

  // Read NZA field
  fscanf (file, "%d %s %s %s %s %s %s", &n, &s1, &s2, &s3, &s4, &s5, &s6);
  Field NZA_1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading NZA_1\n");
	  exit (1);
	}
      else
	{
	  NZA_1.PushData (i, x, y, dydx);
	}
    }

  fclose (file);

  // #################
  // Read second pFile
  // #################
 
  // Check for existence of pfile
  file = OpenFiler (pFile2);

  // Read ne field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ne_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ne_2\n");
	  exit (1);
	}
      else
	{
	  ne_2.PushData (i, x, y, dydx);
	}
    }

  // Read Te field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field Te_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading Te_2\n");
	  exit (1);
	}
      else
	{
	  Te_2.PushData (i, x, y, dydx);
	}
    }

  // Read ni field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ni_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ni_2\n");
	  exit (1);
	}
      else
	{
	  ni_2.PushData (i, x, y, dydx);
	}
    }

  // Read Ti field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field Ti_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading Ti_2\n");
	  exit (1);
	}
      else
	{
	  Ti_2.PushData (i, x, y, dydx);
	}
    }

  // Read nb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field nb_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading nb_2\n");
	  exit (1);
	}
      else
	{
	  nb_2.PushData (i, x, y, dydx);
	}
    }

  // Read pb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field pb_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading pb_2\n");
	  exit (1);
	}
      else
	{
	  pb_2.PushData (i, x, y, dydx);
	}
    }

  // Read ptot field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ptot_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ptot_2\n");
	  exit (1);
	}
      else
	{
	  ptot_2.PushData (i, x, y, dydx);
	}
    }

  // Read omeg field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omeg_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omeg_2\n");
	  exit (1);
	}
      else
	{
	  omeg_2.PushData (i, x, y, dydx);
	}
    }

  // Read omegp field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omegp_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omegp_2\n");
	  exit (1);
	}
      else
	{
	  omegp_2.PushData (i, x, y, dydx);
	}
    }

  // Read omgvb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
 Field omgvb_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgvb_2\n");
	  exit (1);
	}
      else
	{
	  omgvb_2.PushData (i, x, y, dydx);
	}
    }

  // Read omgpp field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omgpp_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgpp_2\n");
	  exit (1);
	}
      else
	{
	  omgpp_2.PushData (i, x, y, dydx);
	}
    }
  
  // Read omgeb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field wE_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgeb_2\n");
	  exit (1);
	}
      else
	{
	  wE_2.PushData (i, x, y, dydx);
	}
    }

  // Read er field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field er_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading er_2\n");
	  exit (1);
	}
      else
	{
	  er_2.PushData (i, x, y, dydx);
	}
    }

  // Read ommvb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ommvb_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ommvb_2\n");
	  exit (1);
	}
      else
	{
	  ommvb_2.PushData (i, x, y, dydx);
	}
    }

  // Read ommpp field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ommpp_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ommpp_2\n");
	  exit (1);
	}
      else
	{
	  ommpp_2.PushData (i, x, y, dydx);
	}
    }

  // Read omevb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omevb_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omevb_2\n");
	  exit (1);
	}
      else
	{
	  omevb_2.PushData (i, x, y, dydx);
	}
    }

  // Read omepp field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omepp_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omepp_2\n");
	  exit (1);
	}
      else
	{
	  omepp_2.PushData (i, x, y, dydx);
	}
    }

  // Read kpol field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field kpol_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading kpol_2\n");
	  exit (1);
	}
      else
	{
	  kpol_2.PushData (i, x, y, dydx);
	}
    }
  
  // Read omghb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omghb_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omghb_2\n");
	  exit (1);
	}
      else
	{
	  omghb_2.PushData (i, x, y, dydx);
	}
    }

  // Read nI field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field nI_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading nz1_2\n");
	  exit (1);
	}
      else
	{
	  nI_2.PushData (i, x, y, dydx);
	}
    }

  // Read vtor1 field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field vtor1_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading vtor1_2\n");
	  exit (1);
	}
      else
	{
	  vtor1_2.PushData (i, x, y, dydx);
	}
    }

  // Read vpol1 field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field vpol1_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading vpol1_2\n");
	  exit (1);
	}
      else
	{
	  vpol1_2.PushData (i, x, y, dydx);
	}
    }

  // Read NZA field
  fscanf (file, "%d %s %s %s %s %s %s", &n, &s1, &s2, &s3, &s4, &s5, &s6);
  Field NZA_2 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading NZA_2\n");
	  exit (1);
	}
      else
	{
	  NZA_2.PushData (i, x, y, dydx);
	}
    }

  fclose (file);

  // ################
  // Read third pFile
  // ################
 
  // Check for existence of pfile
  file = OpenFiler (pFile3);

  // Read ne field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ne_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ne_3\n");
	  exit (1);
	}
      else
	{
	  ne_3.PushData (i, x, y, dydx);
	}
    }

  // Read Te field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field Te_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading Te_3\n");
	  exit (1);
	}
      else
	{
	  Te_3.PushData (i, x, y, dydx);
	}
    }

  // Read ni field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ni_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ni_3\n");
	  exit (1);
	}
      else
	{
	  ni_3.PushData (i, x, y, dydx);
	}
    }

  // Read Ti field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field Ti_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading Ti_3\n");
	  exit (1);
	}
      else
	{
	  Ti_3.PushData (i, x, y, dydx);
	}
    }

  // Read nb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field nb_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading nb_3\n");
	  exit (1);
	}
      else
	{
	  nb_3.PushData (i, x, y, dydx);
	}
    }

  // Read pb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field pb_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading pb_3\n");
	  exit (1);
	}
      else
	{
	  pb_3.PushData (i, x, y, dydx);
	}
    }

  // Read ptot field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ptot_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ptot_3\n");
	  exit (1);
	}
      else
	{
	  ptot_3.PushData (i, x, y, dydx);
	}
    }

  // Read omeg field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omeg_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omeg_3\n");
	  exit (1);
	}
      else
	{
	  omeg_3.PushData (i, x, y, dydx);
	}
    }

  // Read omegp field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omegp_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omegp_3\n");
	  exit (1);
	}
      else
	{
	  omegp_3.PushData (i, x, y, dydx);
	}
    }

  // Read omgvb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omgvb_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgvb_3\n");
	  exit (1);
	}
      else
	{
	  omgvb_3.PushData (i, x, y, dydx);
	}
    }

  // Read omgpp field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omgpp_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgpp_3\n");
	  exit (1);
	}
      else
	{
	  omgpp_3.PushData (i, x, y, dydx);
	}
    }
  
  // Read omgeb field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field wE_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omgeb_3\n");
	  exit (1);
	}
      else
	{
	  wE_3.PushData (i, x, y, dydx);
	}
    }

  // Read er field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field er_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading er_3\n");
	  exit (1);
	}
      else
	{
	  er_3.PushData (i, x, y, dydx);
	}
    }

  // Read ommvb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ommvb_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ommvb_3\n");
	  exit (1);
	}
      else
	{
	  ommvb_3.PushData (i, x, y, dydx);
	}
    }

  // Read ommpp field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field ommpp_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading ommpp_3\n");
	  exit (1);
	}
      else
	{
	  ommpp_3.PushData (i, x, y, dydx);
	}
    }

  // Read omevb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omevb_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omevb_3\n");
	  exit (1);
	}
      else
	{
	  omevb_3.PushData (i, x, y, dydx);
	}
    }

  // Read omepp field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omepp_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omepp_3\n");
	  exit (1);
	}
      else
	{
	  omepp_3.PushData (i, x, y, dydx);
	}
    }

  // Read kpol field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field kpol_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading kpol_3\n");
	  exit (1);
	}
      else
	{
	  kpol_3.PushData (i, x, y, dydx);
	}
    }
  
  // Read omghb field
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field omghb_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading omghb_3\n");
	  exit (1);
	}
      else
	{
	  omghb_3.PushData (i, x, y, dydx);
	}
    }

  // Read nI field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field nI_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading nz1_3\n");
	  exit (1);
	}
      else
	{
	  nI_3.PushData (i, x, y, dydx);
	}
    }

  // Read vtor1 field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field vtor1_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading vtor1_3\n");
	  exit (1);
	}
      else
	{
	  vtor1_3.PushData (i, x, y, dydx);
	}
    }

  // Read vpol1 field 
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  Field vpol1_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading vpol1_3\n");
	  exit (1);
	}
      else
	{
	  vpol1_3.PushData (i, x, y, dydx);
	}
    }

  // Read NZA field
  fscanf (file, "%d %s %s %s %s %s %s", &n, &s1, &s2, &s3, &s4, &s5, &s6);
  Field NZA_3 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("Error reading NZA_3\n");
	  exit (1);
	}
      else
	{
	  NZA_3.PushData (i, x, y, dydx);
	}
    }

  fclose (file);

  // ########################
  // Interpolate profile data
  // ########################
  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);

  Field ne;
  FieldInterpolate (ne_1,    ne_2,    ne_3,    ne,    weight1, weight2, weight3);
  Field Te;
  FieldInterpolate (Te_1,    Te_2,    Te_3,    Te,    weight1, weight2, weight3);
  Field ni;
  FieldInterpolate (ni_1,    ni_2,    ni_3,    ni,    weight1, weight2, weight3);
  Field Ti;
  FieldInterpolate (Ti_1,    Ti_2,    Ti_3,    Ti,    weight1, weight2, weight3);
  Field nb;
  FieldInterpolate (nb_1,    nb_2,    nb_3,    nb,    weight1, weight2, weight3);
  Field pb;
  FieldInterpolate (pb_1,    pb_2,    pb_3,    pb,    weight1, weight2, weight3);
  Field ptot;
  FieldInterpolate (ptot_1,  ptot_2,  ptot_3,  ptot,  weight1, weight2, weight3);
  Field omeg;
  FieldInterpolate (omeg_1,  omeg_2,  omeg_3,  omeg,  weight1, weight2, weight3);
  Field omegp;
  FieldInterpolate (omegp_1, omegp_2, omegp_3, omegp, weight1, weight2, weight3);
  Field omgvb;
  FieldInterpolate (omgvb_1, omgvb_2, omgvb_3, omgvb, weight1, weight2, weight3);
  Field omgpp;
  FieldInterpolate (omgpp_1, omgpp_2, omgpp_3, omgpp, weight1, weight2, weight3);
  Field wE;
  FieldInterpolate (wE_1,    wE_2,    wE_3,    wE,    weight1, weight2, weight3);
  Field er;
  FieldInterpolate (er_1,    er_2,    er_3,    er,    weight1, weight2, weight3);
  Field ommvb;
  FieldInterpolate (ommvb_1, ommvb_2, ommvb_3, ommvb, weight1, weight2, weight3);
  Field ommpp;
  FieldInterpolate (ommpp_1, ommpp_2, ommpp_3, ommpp, weight1, weight2, weight3);
  Field omevb;
  FieldInterpolate (omevb_1, omevb_2, omevb_3, omevb, weight1, weight2, weight3);
  Field omepp;
  FieldInterpolate (omepp_1, omepp_2, omepp_3, omepp, weight1, weight2, weight3);
  Field kpol;
  FieldInterpolate (kpol_1,  kpol_2,  kpol_3,  kpol,  weight1, weight2, weight3);
  Field omghb;
  FieldInterpolate (omghb_1, omghb_2, omghb_3, omghb, weight1, weight2, weight3);
  Field nI;
  FieldInterpolate (nI_1,    nI_2,    nI_3,    nI,    weight1, weight2, weight3);
  Field vtor1;
  FieldInterpolate (vtor1_1, vtor1_2, vtor1_3, vtor1, weight1, weight2, weight3);
  Field vpol1;
  FieldInterpolate (vpol1_1, vpol1_2, vpol1_3, vpol1, weight1, weight2, weight3);
  Field NZA;
  FieldInterpolate (NZA_1,   NZA_2,   NZA_3,   NZA,   weight1, weight2, weight3);

  // ########################
  // Write interpolated pFile
  // ########################
  file = OpenFilew (pFile);

  fprintf (file, "%3d %s %s %s", ne.GetN (), "PsiN", " ne", " dne/dPsiN\n");
  for (int i = 0; i < ne.GetN (); i++)
    {
      ne.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Te.GetN (), "PsiN", " Te", " dTe/dPsiN\n");
  for (int i = 0; i < Te.GetN (); i++)
    {
      Te.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ni.GetN (), "PsiN", " ni", " dni/dPsiN\n");
  for (int i = 0; i < ni.GetN (); i++)
    {
      ni.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Ti.GetN (), "PsiN", " Ti", " dTi/dPsiN\n");
  for (int i = 0; i < Ti.GetN (); i++)
    {
      Ti.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nb.GetN (), "PsiN", " nb", " dnb/dPsiN\n");
  for (int i = 0; i < nb.GetN (); i++)
    {
      nb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", pb.GetN (), "PsiN", " pb", " dpb/dPsiN\n");
  for (int i = 0; i < pb.GetN (); i++)
    {
      pb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ptot.GetN (), "PsiN", " ptot", " dptot/dPsiN\n");
  for (int i = 0; i < ptot.GetN (); i++)
    {
      ptot.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omeg.GetN (), "PsiN", " omeg", " domeg/dPsiN\n");
  for (int i = 0; i < omeg.GetN (); i++)
    {
      omeg.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omegp.GetN (), "PsiN", " omegp", " domegp/dPsiN\n");
  for (int i = 0; i < omegp.GetN (); i++)
    {
      omegp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omgvb.GetN (), "PsiN", " omgvb", " domgvb/dPsiN\n");
  for (int i = 0; i < omgvb.GetN (); i++)
    {
      omgvb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omgpp.GetN (), "PsiN", " omgpp", " domgpp/dPsiN\n");
  for (int i = 0; i < omgpp.GetN (); i++)
    {
      omgpp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wE.GetN (), "PsiN", " wE", " dwE/dPsiN\n");
  for (int i = 0; i < wE.GetN (); i++)
    {
      wE.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", er.GetN (), "PsiN", " er", " der/dPsiN\n");
  for (int i = 0; i < er.GetN (); i++)
    {
      er.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ommvb.GetN (), "PsiN", " ommvb", " dommvb/dPsiN\n");
  for (int i = 0; i < ommvb.GetN (); i++)
    {
      ommvb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ommpp.GetN (), "PsiN", " ommpp", " dommpp/dPsiN\n");
  for (int i = 0; i < ommpp.GetN (); i++)
    {
      ommpp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omevb.GetN (), "PsiN", " omevb", " domevb/dPsiN\n");
  for (int i = 0; i < omevb.GetN (); i++)
    {
      omevb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omepp.GetN (), "PsiN", " omepp", " domepp/dPsiN\n");
  for (int i = 0; i < omepp.GetN (); i++)
    {
      omepp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", kpol.GetN (), "PsiN", " kpol", " dkpol/dPsiN\n");
  for (int i = 0; i < kpol.GetN (); i++)
    {
      kpol.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", omghb.GetN (), "PsiN", " omghb", " domghb/dPsiN\n");
  for (int i = 0; i < omghb.GetN (); i++)
    {
      omghb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nI.GetN (), "PsiN", " nI", " dnI/dPsiN\n");
  for (int i = 0; i < nI.GetN (); i++)
    {
      nI.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", vtor1.GetN (), "PsiN", " vtor1", " dvtor1/dPsiN\n");
  for (int i = 0; i < vtor1.GetN (); i++)
    {
      vtor1.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", vpol1.GetN (), "PsiN", " vpol1", " dvpol1/dPsiN\n");
  for (int i = 0; i < vpol1.GetN (); i++)
    {
      vpol1.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", NZA.GetN (), "N", " Z", " A of ION SPECIES\n");
  for (int i = 0; i < NZA.GetN (); i++)
    {
      NZA.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
    
  fclose (file);

  printf ("pFile Interpolation:\n");
  printf ("%s %11.4e\n", pFile1, weight1);
  printf ("%s %11.4e\n", pFile2, weight2);
  printf ("%s %11.4e\n", pFile3, weight3);
}

void Neoclassical::pFileInterp (vector<string> pFileName, vector<double> pFileTime, int pFileNumber, double time)
{
  int    index;
  double _time;

  if (time < pFileTime[0])
    {
      index = 0;
      _time = pFileTime[0];
    }
  else if (time >= pFileTime[pFileNumber-1])
    {
      index = pFileNumber - 3;
      _time = pFileTime[pFileNumber-1];
    }
  else
    {
      for (int i = 0; i < pFileNumber-1; i++)
	if (time >= pFileTime[i] && time < pFileTime[i+1])
	  {
	    index = i;
	    _time = time;

	    if (index > pFileNumber-3)
	      index = pFileNumber - 3;
	  }
    }

  char* pFile = "pFile";
  char* file1 = (char*) pFileName[index  ].c_str();
  char* file2 = (char*) pFileName[index+1].c_str();
  char* file3 = (char*) pFileName[index+2].c_str();
  
  pFileInterpolate (file1, pFileTime[index], file2, pFileTime[index+1], file3, pFileTime[index+2], pFile, _time);
}
