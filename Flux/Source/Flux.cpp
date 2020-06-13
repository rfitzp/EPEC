// Flux.cpp

// PROGRAM ORGANIZATION:
//
//       Flux:: Flux          ()       
// void  Flux:: Solve         (int _INTP, int _NTOR, int _MMIN, int _MMAX, double _TIME)
// void  Flux:: SetParameters (int _INTG, int _NTOR, int _MMIN, int _MMAX, double _TIME)
// FILE* Flux:: OpenFilew     (char* filename)
// FILE* Flux:: OpenFiler     (char* filename)
// FILE* Flux:: OpenFilea     (char* filename)

#include "Flux.h"

// ###########
// Constructor
// ###########
Flux::Flux ()
{
}

// #########################
// Function to solve problem
// #########################
void Flux::Solve (int _INTP, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM)
{
  // Set global parameters
  SetParameters (_INTP, _NTOR, _MMIN, _MMAX, _TIME, _PSILIM);

  // Input gFile data and output Stage1 data.
  // Stage1 data output to directory /Stage1.
  Stage1 ();
  
  // Input Stage1 data and output Stage2 data.
  // Stage2 data output to directory /Stage2.
  // Data passed to other programs output to fFile.
  Stage2 ();
}

// #################################
// Function to set global parameters
// #################################
void Flux::SetParameters (int _INTG, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM)
{
  // Set default values of input parameters
  INTG   = 0;
  NPSI   = 256;
  NTHETA = 512;
  NNC    = 10;
  QFLG   = 0;
  Q95    = 2.5;
  NTOR   = 2;
  MMIN   = 2;
  MMAX   = 20;
  PSILIM = 0.997;
  TIME   = 0.;
  
  H0     = 1.e-6;
  ACC    = 1.e-14;
  ETA    = 1.e-8;
  DR     = 1.e-2;
 
  // Read namelist
  NameListRead (&INTG, &NPSI, &NTHETA, &NNC, &NTOR, &QFLG, &Q95, &H0, &ACC, &ETA, &DR, &MMIN, &MMAX, &PSILIM, &TIME);

  // Override namelist values with command line options
  if (_NTOR > 0)
    NTOR = _NTOR;
  if (_MMIN > 0)
    MMIN = _MMIN;
  if (_MMAX > 0)
    MMAX = _MMAX;
  if (_TIME > 0.)
    TIME = _TIME;
  if (_INTG > -1)
    INTG = _INTG;
  if (_PSILIM > 0.)
    PSILIM = _PSILIM;

  // Output calculation parameters
  printf ("Input Parameters (from Inputs/Flux.in and command line options):\n");
  printf ("NPSI = %4d         NTHETA = %4d         NNC  = %3d\n",
	  NPSI, NTHETA, NNC);
  printf ("Q95  = %11.4e  QFLG   = %2d\n",
	  Q95, QFLG, NTOR);
  printf ("NTOR = %2d           MMIN   = %2d           MMAX =  %2d          PSILIM = %11.4e  TIME = %11.4e  INTG = %2d\n",
	  NTOR, MMIN, MMAX, PSILIM, TIME, INTG);
  printf ("H0   = %11.4e  ACC    = %11.4e  ETA  = %11.4e      DR = %11.4e\n",
	  H0, ACC, ETA, DR);

  FILE* namelist = OpenFilew ((char*) "Inputs/InputParameters.txt");
  fprintf (namelist, "Input Parameters (from Inputs/Flux.in and command line options):\n");
  fprintf (namelist, "NPSI = %4d         NTHETA = %4d         NNC  = %3d\n",
	  NPSI, NTHETA, NNC);
  fprintf (namelist, "Q95  = %11.4e  QFLG   = %2d\n",
	  Q95, QFLG, NTOR);
  fprintf (namelist, "NTOR = %2d           MMIN   = %2d           MMAX =  %2d          PSILIM = %11.4e  TIME = %11.4e  INTG = %2d\n",
	  NTOR, MMIN, MMAX, PSILIM, TIME, INTG);
  fprintf (namelist, "H0   = %11.4e  ACC    = %11.4e  ETA  = %11.4e      DR = %11.4e\n",
	  H0, ACC, ETA, DR);
  fclose (namelist);
  
  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "Input Parameters (from Inputs/Flux.in and command line options):\n");
  fprintf (monitor, "NPSI = %4d         NTHETA = %4d         NNC  = %3d\n",
	  NPSI, NTHETA, NNC);
  fprintf (monitor, "Q95  = %11.4e  QFLG   = %2d\n",
	  Q95, QFLG, NTOR);
  fprintf (monitor, "NTOR = %2d           MMIN   = %2d           MMAX =  %2d          PSILIM = %11.4e  TIME = %11.4e  INTG = %2d\n",
	  NTOR, MMIN, MMAX, PSILIM, TIME, INTG);
  fprintf (monitor, "H0   = %11.4e  ACC    = %11.4e  ETA  = %11.4e      DR = %11.4e\n",
	  H0, ACC, ETA, DR);
  fclose (monitor);

  // Sanity check
  if (NPSI < 1)
    {
      printf ("FLUX::SetParameters: Error - NPSI must be positive\n");
      exit (1);
    } 
  if (NTHETA < 1)
    {
      printf ("FLUX::SetParameters: Error NTHETA must be positive\n");
      exit (1);
    }
  if (NTHETA%2 == 0)
    {
      printf ("FLUX::SetParameters: Error - NTHETA must be odd\n");
      exit (1);
    }
  if (NNC < 1)
    {
      printf ("FLUX::SetParameters: Error - NNC must be positive\n");
      exit (1);
    } 
  if (NTOR < 1)
    {
      printf ("FLUX::SetParameters: Error - NTOR must be positive\n");
      exit (1);
    }
  if (MMIN < 1)
    {
      printf ("FLUX::SetParameters: Error - MMIN must be positive\n");
      exit (1);
    }
  if (MMAX < 1)
    {
      printf ("FLUX::SetParameters: Error - MMAX must be positive\n");
      exit (1);
    }
  if (MMAX <= MMIN)
    {
      printf ("FLUX::SetParameters: Error - MMAX must be greater than MMIN\n");
      exit (1);
    }
  if (PSILIM <= 0. || PSILIM > 1.)
    {
      printf ("FLUX::SetParameters: Error - PSILIM must lie between 0 and 1\n");
      exit (1);
    }
  if (QFLG && Q95 <= 0.)
    {
      printf ("FLUX::SetParameters: Error - Q95 must be positive\n");
      exit (1);
    }
  
  FILE* file = OpenFilew ((char *) "Outputs/Stage2/NpsiNtor.txt");
  fprintf (file, "%d %d\n", NPSI, NTOR);
  fclose (file);
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Flux::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("FLUX::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for reading
// ##########################################
FILE* Flux::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("FLUX::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ############################################
// Function to open existing file for appending
// ############################################
FILE* Flux::OpenFilea (char* filename)
{
  FILE* file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("FLUX::OpenFilea: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

