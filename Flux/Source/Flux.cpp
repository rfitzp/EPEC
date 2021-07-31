// Flux.cpp

// PROGRAM ORGANIZATION:
//
//       Flux:: Flux          ()       
// void  Flux:: Solve         (int _INTP, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM, double _PSIPED, double _PSIRAT)
// void  Flux:: SetParameters (int _INTG, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM, double _PSIPED, double _PSIRAT)
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
void Flux::Solve (int _INTP, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM, double _PSIPED, double _PSIRAT)
{
  // Set global parameters
  SetParameters (_INTP, _NTOR, _MMIN, _MMAX, _TIME, _PSILIM, _PSIPED, _PSIRAT);

  // Input gFile data and output Stage1 data.
  // Stage1 data output to directory Outputs/Stage1.
  Stage1 ();
  
  // Input Stage1 data and output Stage2 data.
  // Stage2 data output to directory Outputs/Stage2.
  // Data passed to other programs output to Outputs/fFile.
  Stage2 ();
}

// #################################
// Function to set global parameters
// #################################
void Flux::SetParameters (int _INTG, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM, double _PSIPED, double _PSIRAT)
{
  // Set default values of input parameters
  NTOR    = 2;
  MMIN    = 2;
  MMAX    = 20;

  PSILIM  = 0.998;
  PSIRAT  = 0.995;
  PSIPED  = 0.950;

  INTG    = 0;
  TIME    = 0.;
  
  NPSI    = 256;
  PACK    = 1.;
  NTHETA  = 512;
  NNC     = 10;
  NSMOOTH = 100;

  H0      = 1.e-6;
  ACC     = 1.e-14;
  ETA     = 1.e-8;
  
  // Read namelist file Inputs/Flux.nml
  NameListRead (&INTG, &NPSI, &PACK, &NTHETA, &NNC, &NTOR, &H0, &ACC, &ETA, &MMIN, &MMAX, &PSILIM, &TIME, &PSIPED, &NSMOOTH, &PSIRAT);

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
  if (_PSIPED > 0.)
    PSIPED = _PSIPED;
  if (_PSIRAT > 0.)
    PSIRAT = _PSIRAT;

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
  if (PSIRAT < 0. || PSIRAT > PSILIM)
    {
      printf ("FLUX::SetParameters: Error - PSIRAT cannot be negative or exceed PSILIM\n");
      exit (1);
    }
  if (H0 <= 0.)
    {
      printf ("FLUX::SetParameters: Error - H0 must be positive\n");
      exit (1);
    }
  if (ACC <= 0.)
    {
      printf ("FLUX::SetParameters: Error - ACC must be positive\n");
      exit (1);
    }
  if (ETA <= 0.)
    {
      printf ("FLUX::SetParameters: Error - ETA must be positive\n");
      exit (1);
    }
  if (PSIPED <= 0. || PSIPED > 1.)
    {
      printf ("FLUX::SetParameters: Error - PSIPED must lie betweeen 0 and 1\n");
      exit (1);
    }
 
   // Output PSIPED
   FILE* Pfile = OpenFilew ((char*) "Outputs/Stage1/Psilim.txt");
   fprintf (Pfile, "%16.9e %16.9e %16.9e\n", PSILIM, PSIPED, PSIRAT);
   fclose (Pfile);
   
   // Output calculation parameters
   printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
   printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
   printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
   printf ("Input Parameters (from Inputs/Flux.nml and command line options):\n");
   printf ("NPSI = %4d         NTHETA = %4d         NNC  = %3d          PACK   = %11.4e\n",
	   NPSI, NTHETA, NNC, PACK);
   printf ("NTOR = %2d           MMIN   = %2d           MMAX =  %2d          PSILIM = %11.4e  TIME   = %11.4e  INTG = %2d  PSIPED = %11.4e  NSMOOTH = %3d  PSIRAT = %11.4e\n",
	   NTOR, MMIN, MMAX, PSILIM, TIME, INTG, PSIPED, NSMOOTH, PSIRAT);
   printf ("H0   = %11.4e  ACC    = %11.4e  ETA  = %11.4e\n",
	   H0, ACC, ETA);
   
   FILE* namelist = OpenFilew ((char*) "Inputs/InputParameters.txt");
   fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
   fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
   fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
   fprintf (namelist, "Input Parameters (from Inputs/Flux.nml and command line options):\n");
   fprintf (namelist, "NPSI = %4d         NTHETA = %4d         NNC  = %3d          PACK   = %11.4e\n",
	    NPSI, NTHETA, NNC, PACK);
   fprintf (namelist, "NTOR = %2d           MMIN   = %2d           MMAX =  %2d          PSILIM = %11.4e  TIME   = %11.4e  INTG = %2d  PSIPED = %11.4e  NSMOOTH = %11.4e  PSIRAT = %11.4e\n",
	    NTOR, MMIN, MMAX, PSILIM, TIME, INTG, PSIPED, NSMOOTH, PSIRAT);
   fprintf (namelist, "H0   = %11.4e  ACC    = %11.4e  ETA  = %11.4e\n",
	    H0, ACC, ETA);
   fclose (namelist);
   
   FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
   fprintf (monitor, "Git Hash     = "); fprintf (monitor, GIT_HASH);     fprintf (monitor, "\n");
   fprintf (monitor, "Compile time = "); fprintf (monitor, COMPILE_TIME); fprintf (monitor, "\n");
   fprintf (monitor, "Git Branch   = "); fprintf (monitor, GIT_BRANCH);   fprintf (monitor, "\n\n");
   fprintf (monitor, "Input Parameters (from Inputs/Flux.nml and command line options):\n");
   fprintf (monitor, "NPSI = %4d         NTHETA = %4d         NNC  = %3d          PACK   = %11.4e\n",
	    NPSI, NTHETA, NNC, PACK);
   fprintf (monitor, "NTOR = %2d           MMIN   = %2d           MMAX =  %2d          PSILIM = %11.4e  TIME   = %11.4e  INTG = %2d  PSIPED = %11.4e  NSMOOTH = %11.4e  PSIRAT = %11.4e\n",
	    NTOR, MMIN, MMAX, PSILIM, TIME, INTG, PSIPED, NSMOOTH, PSIRAT);
   fprintf (monitor, "H0   = %11.4e  ACC    = %11.4e  ETA  = %11.4e\n",
	    H0, ACC, ETA);
   fclose (monitor);
   
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

// #################################
// Function to call operating system
// #################################
void Flux::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("FLUX: Operating system call error executing %s\n", command);
      exit (1);
    }
}
