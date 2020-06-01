// pFileRead.cpp

#include "Neoclassical.h"
#include "Field.h"

// ###############################
// Function to read standard pfile
// ###############################
void Neoclassical::pFileRead ()
{
  int    n;
  double x, y, dydx;
  char   s1[30], s2[30], s3[30], s4[30], s5[30], s6[30];
  
  // Check for existence of pfile
  FILE* file = OpenFiler ((char*) "pFile");
  if (file == NULL) 
    {
      printf ("NEOCLASSICAL::pFileRread: Error opening pFile\n");
      exit (1);
    }

  // Read ne field (assumed units - 10^20 m^-3)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading ne    from pFile - n = %4d:\n", n);
  ne.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading ne\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e20;
	  dydx *= 1.e20;
	  ne.PushData (i, x, y, dydx);
	}
    }

  // Read Te field (assumed units - keV)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading Te    from pFile - n = %4d:\n", n);
  Te.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading Te\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e3 * e;
	  dydx *= 1.e3 * e;
	  Te.PushData (i, x, y, dydx);
	}
    }

  // Read ni field (assumed units - 10^20 m^-3)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading ni    from pFile - n = %4d:\n", n);
  ni.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading ni\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e20;
	  dydx *= 1.e20;
	  ni.PushData (i, x, y, dydx);
	}
    }

  // Read Ti field (assumed units - keV)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading Ti    from pFile - n = %4d:\n", n);
  Ti.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading Ti\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e3 * e;
	  dydx *= 1.e3 * e;
	  Ti.PushData (i, x, y, dydx);
	}
    }

  // Read nb field (assumed units - 10^20 m^-3)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading nb    from pFile - n = %4d:\n", n);
  nb.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading nb\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e20;
	  dydx *= 1.e20;
	  nb.PushData (i, x, y, dydx);
	}
    }

  // Read pb field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading pb    from pFile - n = %4d:\n", n);
  Field pb (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading pb\n");
	  exit (1);
	}
      else
	{
	  pb.PushData (i, x, y, dydx);
	}
    }

  // Read ptot field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading ptot  from pFile - n = %4d:\n", n);
  Field ptot (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading ptot\n");
	  exit (1);
	}
      else
	{
	  ptot.PushData (i, x, y, dydx);
	}
    }

  // Read omeg field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omeg  from pFile - n = %4d:\n", n);
  Field omeg (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omeg\n");
	  exit (1);
	}
      else
	{
	  omeg.PushData (i, x, y, dydx);
	}
    }

  // Read omegp field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omegp from pFile - n = %4d:\n", n);
  Field omegp (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omegp\n");
	  exit (1);
	}
      else
	{
	  omegp.PushData (i, x, y, dydx);
	}
    }

  // Read omgvb field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omgvb from pFile - n = %4d:\n", n);
  Field omgvb (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omgvb\n");
	  exit (1);
	}
      else
	{
	  omgvb.PushData (i, x, y, dydx);
	}
    }

  // Read omgpp field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omgpp from pFile - n = %4d:\n", n);
  Field omgpp (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omgpp\n");
	  exit (1);
	}
      else
	{
	  omgpp.PushData (i, x, y, dydx);
	}
    }
  
  // Read omgeb field (assumed units krad/s)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omgeb from pFile - n = %4d:\n", n);
  wE.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omgeb\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e3;
	  dydx *= 1.e3;
	  wE.PushData (i, x, y, dydx);
	}
    }

  // Read er field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading er    from pFile - n = %4d:\n", n);
  Field er (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading er\n");
	  exit (1);
	}
      else
	{
	  er.PushData (i, x, y, dydx);
	}
    }

  // Read ommvb field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading ommvb from pFile - n = %4d:\n", n);
  Field ommvb (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading ommvb\n");
	  exit (1);
	}
      else
	{
	  ommvb.PushData (i, x, y, dydx);
	}
    }

  // Read ommpp field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading ommpp from pFile - n = %4d:\n", n);
  Field ommpp (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading ommpp\n");
	  exit (1);
	}
      else
	{
	  ommpp.PushData (i, x, y, dydx);
	}
    }

  // Read omevb field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omevb from pFile - n = %4d:\n", n);
  Field omevb (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omevb\n");
	  exit (1);
	}
      else
	{
	  omevb.PushData (i, x, y, dydx);
	}
    }

  // Read omepp field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omepp from pFile - n = %4d:\n", n);
  Field omepp (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omepp\n");
	  exit (1);
	}
      else
	{
	  omepp.PushData (i, x, y, dydx);
	}
    }

  // Read kpol field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading kpol  from pFile - n = %4d:\n", n);
  Field kpol (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading kpol\n");
	  exit (1);
	}
      else
	{
	  kpol.PushData (i, x, y, dydx);
	}
    }
  
  // Read omghb field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading omghb from pFile - n = %4d:\n", n);
  Field omghb (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading omghb\n");
	  exit (1);
	}
      else
	{
	  omghb.PushData (i, x, y, dydx);
	}
    }

  // Read nI1 field (assumed units 10^20 m^-3)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading nI1   from pFile - n = %4d:\n", n);
  nI.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading nI1\n");
	  exit (1);
	}
      else
	{
	  y    *= 1.e20;
	  dydx *= 1.e20;
	  nI.PushData (i, x, y, dydx);
	}
    }

  // Read vtor1 field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading vtor1 from pFile - n = %4d:\n", n);
  Field vtor1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading vtor1\n");
	  exit (1);
	}
      else
	{
	  vtor1.PushData (i, x, y, dydx);
	}
    }

  // Read vpol1 field (not used)
  fscanf (file, "%d %s %s %s", &n, &s1, &s2, &s3);
  printf ("Reading vpol1 from pFile - n = %4d:\n", n);
  Field vpol1 (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading vpol1\n");
	  exit (1);
	}
      else
	{
	  vpol1.PushData (i, x, y, dydx);
	}
    }

  // Read NZA field
  fscanf (file, "%d %s %s %s %s %s %s", &n, &s1, &s2, &s3, &s4, &s5, &s6);
  printf ("Reading NZA   from pFile - n = %4d:\n", n);
  NZA.resize (n);
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading NZA\n");
	  exit (1);
	}
      else
	{
	  NZA.PushData (i, x, y, dydx);
	}
    }

  fclose (file);
}


