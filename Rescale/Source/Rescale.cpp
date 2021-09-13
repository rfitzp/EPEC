// Rescale.cpp

#include "Rescale.h"

// ###########
// Constructor
// ###########
Rescale::Rescale ()
{
}

// ###############################
// Function to rescale equilibrium
// ###############################
void Rescale::RescaleEquilibrium ()
{
  // Rescale gFile
  gFileRescale (&q95_old, &q95_new, &a1);

  printf ("\nq95_old = %11.4e  q95_new = %11.4e  rescale factor = %11.4e\n\n", q95_old, q95_new, a1);
  
  // Rescale pFile
  pFileRescale ();
  printf ("\n");
}

// #########################
// Function to rescale pFile
// #########################
void Rescale::pFileRescale ()
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescale: Error opening input pFile\n");
      exit (1);
    }

  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescale: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("NEOCLASSICAL::pFileRead: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double A1 = a1;
      if (strstr (s, "ION") != NULL)
	A1 = 1.;
      else if (strstr (s, "pb") != NULL)
	A1 = a1*a1;
      else if (strstr (s, "ptot") != NULL)
	A1 = a1*a1;

      field.Rescale (A1);

      printf ("pFile field: %-34s:  rescale factor = %11.4e\n", s, A1);
    }
  while (1);
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Rescale::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("RESCALE::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for reading
// ##########################################
FILE* Rescale::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("RESCALE::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ############################################
// Function to open existing file for appending
// ############################################
FILE* Rescale::OpenFilea (char* filename)
{
  FILE* file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("RESCALE::OpenFilea: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to call operating system
// #################################
void Rescale::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("RESCALE: Operating system call error executing %s\n", command);
      exit (1);
    }
}

