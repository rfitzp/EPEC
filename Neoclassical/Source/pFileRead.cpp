// pFileRead.cpp

#include "Neoclassical.h"

// ###############################
// Function to read standard pfile
// ###############################
void Neoclassical::pFileRead ()
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  
  // Check for existence of pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("NEOCLASSICAL::pFileRread: Error opening pFile\n");
      exit (1);
    }

  // Read data from pFile
  int ne_flag = 0; int Te_flag = 0; int ni_flag = 0; int Ti_flag = 0;
  int nb_flag = 0; int wE_flag = 0; int nI_flag = 0; int NZ_flag = 0;
  int wt_flag = 0; 
  
  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileRead: Error reading pFile\n");
	  exit (1);
	}
      
      if (strstr (s, "ne(") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag = 1;
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
	}
      else if (strstr (s, "te(") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag = 1;
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
	}
      else if (strstr (s, "ni(") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag = 1;
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
	}
      else if (strstr (s, "ti(") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag = 1;
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
	}
      else if (strstr (s, "nb(") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag = 1;
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
	}
      else if (strstr (s, "omgeb(") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag = 1;
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
	}
       else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag = 1;
	  printf ("Reading  omeg from pFile - n = %4d:\n", n);
	  wt.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileRead: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  y    *= 1.e3;
		  dydx *= 1.e3;
		  wt.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1(") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag = 1;
	  printf ("Reading nz1   from pFile - n = %4d:\n", n);
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
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag = 1;
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
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileRead: Error reading unused field\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag * Te_flag * ni_flag * Ti_flag * nb_flag * wE_flag * wt_flag * nI_flag * NZ_flag == 0)
    {
      printf ("NEOCLASSICAL::pFileRead: Missing field in pFile\n");
      exit (1);
    }
}


