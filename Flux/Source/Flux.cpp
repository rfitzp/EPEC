// Flux.cpp

// PROGRAM ORGANIZATION:
//
//       Flux:: Flux          ()       
// void  Flux:: Solve         ()
// void  Flux:: SetParameters ()
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
void Flux::Solve ()
{
  // Set global parameters
  SetParameters ();

  // Input gFile data and output Stage1 data.
  // Stage1 data output to directory Outputs/Stage1.
  Stage1 ();
  
  // Input Stage1 data and output Stage2 data.
  // Stage2 data output to directory Outputs/Stage2.
  // Data passed to other programs output to Outputs/fFile.
  if (STAGE2) Stage2 ();
}

// #################################
// Function to set global parameters
// #################################
void Flux::SetParameters ()
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

  RW      = 1.2;
  
  NPSI    = 256;
  PACK    = 1.;
  NTHETA  = 512;
  NNC     = 10;
  NSMOOTH = 100;
  NEOANG  = 0;

  STAGE2  = 1;

  H0      = 1.e-6;
  ACC     = 1.e-14;
  ETA     = 1.e-8;

  EPS     = 1.e-4;
  DELTA   = 1.e-7;
  
  // Read namelist file Inputs/Flux.nml
  NameListRead (&INTG, &NPSI, &PACK, &NTHETA, &NNC, &NTOR, &H0, &ACC, &ETA, &MMIN, &MMAX, &PSILIM, &TIME, &PSIPED, &NSMOOTH, &PSIRAT, &NEOANG, &RW, &STAGE2);

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
  if (RW <= 1.)
     {
      printf ("FLUX::SetParameters: Error - RW must be greater than unity\n");
      exit (1);
    }
 
   // Output calculation parameters
   printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
   printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
   printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
   printf ("Input Parameters (from Inputs/Flux.nml):\n");
   printf ("NPSI = %4d  NTHETA = %4d  NNC = %3d  PACK = %10.3e  NEOANG = %2d  STAGE2 = %2d\n",
	   NPSI, NTHETA, NNC, PACK, NEOANG, STAGE2);
   printf ("NTOR = %2d  MMIN = %2d  MMAX = %2d  PSILIM = %10.3e  PSIRAT = %10.3e  PSIPED = %10.3e  TIME = %10.3e  INTG = %2d  RW = %10.3e\n",
	   NTOR, MMIN, MMAX, PSILIM, PSIRAT, PSIPED, TIME, INTG, RW);
   printf ("H0 = %10.3e  ACC = %10.3e  ETA = %10.3e  NSMOOTH = %3d\n",
	   H0, ACC, ETA, NSMOOTH);
   
   FILE* namelist = OpenFilew ((char*) "Outputs/InputParameters.txt");
   fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
   fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
   fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
   fprintf (namelist, "Input Parameters (from Inputs/Flux.nml):\n");
   fprintf (namelist, "NPSI = %4d  NTHETA = %4d  NNC = %3d  PACK = %10.3e  NEOANG = %2d  STAGE2 = %2d\n",
	    NPSI, NTHETA, NNC, PACK, NEOANG, STAGE2);
   fprintf (namelist, "NTOR = %2d  MMIN = %2d  MMAX = %2d  PSILIM = %10.3e  PSIRAT = %10.3e  PSIPED = %10.3e  TIME = %10.3e  INTG = %2d  RW = %10.3e\n",
	    NTOR, MMIN, MMAX, PSILIM, PSIRAT, PSIPED, TIME, INTG, RW);
   fprintf (namelist, "H0 = %10.3e  ACC = %10.3e  ETA = %10.3e  NSMOOTH = %3d\n",
	    H0, ACC, ETA, NSMOOTH);
   fclose (namelist);
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
