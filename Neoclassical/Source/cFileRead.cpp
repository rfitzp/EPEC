// cFileRead.cpp

#include "Neoclassical.h"

// ###############################
// Function to read standard cfile
// ###############################
void Neoclassical::cFileRead ()
{
  int    n;
  double x, y, dydx;
 
  // Check for existence of cfile
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

  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf", &x, &y) != 2)
	{
	  printf ("NEOCLASSICAL::cFileRead: Error reading cFile (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip.PushData (i, x, y, dydx);
	}
    }
  
  fclose (file);
}


