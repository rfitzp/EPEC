// cFileRead.cpp

#include "Neoclassical.h"

// ###############################
// Function to read standard cfile
// ###############################
void Neoclassical::cFileRead ()
{
  int    n;
  double x, y1, y2, y3, y4, dydx;
 
  // Check for existence of cFile
  FILE* file = OpenFiler ((char*) "Inputs/cFile");
  if (file == NULL) 
    {
      printf ("NEOCLASSICAL::cFileRead: Error opening cFile\n");
      exit (1);
    }

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileRead: Error reading cFile (1)\n");
      exit (1);
    }

  Chip.resize (n);
  Chie.resize (n);
  Chin.resize (n);
  Chii.resize (n);

  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileRead: Error reading cFile (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip.PushData (i, x, y1, dydx);
	  Chie.PushData (i, x, y2, dydx);
	  Chin.PushData (i, x, y3, dydx);
	  Chii.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);
}


