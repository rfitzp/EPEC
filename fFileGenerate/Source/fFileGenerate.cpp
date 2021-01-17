// fFileGenerate.cpp

// ##############################################################
// Program to generate series of interpolated fFiles and nFiles.

// Uses FLUX/NEOCLASSICAL codes.

// .........
// Versions:
// .........

// 1.0 - Initial version

// ##############################################################

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

#define MAXCOMMANDLINELENGTH 500

extern "C" void NameListRead (int* FLUX_NTOR, int* FLUX_MMIN, int* FLUX_MMAX,
			      int* NEO_IMPURITY, int* NEO_NEUTRAL, int* NEO_FREQ, int* NEO_NTYPE, double* NEO_NN, double* NEO_LN, double* NEO_YN,
			      double* TSTART, double* TEND, double* DT); 

void fFileGenerate ();

// #############
// Main function
// #############

int main ()
{
  fFileGenerate ();

  return 0;
}

// #################################
// Function to call operating system
// #################################
double CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("Operating system call error executing %s\n", command);
      return 1;
    }

  return 0;
}

// ########################
// Function to perform scan
// ########################
void fFileGenerate ()
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n#####################\n");
  printf ("Program fFileGenerate\n");
  printf ("#####################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // .............
  // Read namelist
  // .............
  int    FLUX_NTOR;      // Set toroidal mode number in FLUX
  int    FLUX_MMIN;      // Set minimum poloidal mode number in FLUX
  int	 FLUX_MMAX;      // Set maximum poloidal mode number in FLUX
  int    FLUX_INTG = 1;  // Interpolate gFiles in FLUX
  int    NEO_INTP  = 1;  // Interpolate pFiles in NEOCLASSICAL
  int	 NEO_INTF  = 1;  // Interpolate fFiles in NEOCLASSICAL
  int	 NEO_IMPURITY;   // If != 0 then include single ion impurity species in NEOCLASSICAL
  int	 NEO_NEUTRAL;    // If != 0  then include majority ion neutrals in NEOCLASSICAL
  int	 NEO_FREQ;       // Set natural frequency type (see NEOCLASSICAL/Inputs/Neoclassical.nml)
  int    NEO_NTYPE;      // If 0/1 then neutral density distribution exponential/Lorentzian in NEOCLASSICAL
  double NEO_NN;         // Set neutral density at boundary in NEOCLASSICAL
  double NEO_LN;         // Set neutral density scalelength in NEOCLASSICAL
  double NEO_YN;         // Set neutral poloidal asymmetry parameter in NEOCLASSICAL
  double TSTART;         // Scan start time (ms)
  double TEND;           // Scan end time (ms)
  double DT;             // Scan time step (ms)

  NameListRead (&FLUX_NTOR, &FLUX_MMIN, &FLUX_MMAX,
		&NEO_IMPURITY, &NEO_NEUTRAL, &NEO_FREQ, &NEO_NTYPE, &NEO_NN, &NEO_LN, &NEO_YN,
		&TSTART, &TEND, &DT);

  printf ("Reading Inputs/fFile.nml:\n");
  printf ("FLUX_NTOR = %2d  FLUX_MMIN = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  printf ("NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  printf ("TSTART = %11.4e  TEND = %11.4e  DT = %11.4e\n",
	  TSTART, TEND, DT);

  // ............
  // Sanity check
  // ............
  if (TEND < TSTART)
    {
      printf ("Error - TEND must be greater than TSTART\n");
      exit (1);
    }
  if (DT < 0.)
    {
      printf ("Error - DT must be positive\n");
      exit (1);
    }

  // .........................
  // Initialize symbolic links
  // .........................
  CallSystem ("cd ../Flux/Outputs/fFiles; rm -rf *");
  CallSystem ("cd ../Neoclassical/Outputs/nFiles; rm -rf *");

  // .................
  // Get date and time
  // .................
  time_t     rawtime;
  struct tm* timeinfo;
  time                 (&rawtime);
  timeinfo = localtime (&rawtime);
  
  printf ("%s\n", asctime (timeinfo));
  printf ("\n#####################\n");
  printf ("Program fFileGenerate\n");
  printf ("#####################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  
  printf ("Input parameters (from Inputs/fFile.nml):\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("FLUX_NTOR = %2d  FLUX_MMIN = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  printf ("NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  printf ("TSTART = %11.4e  TEND = %11.4e  DT = %11.4e\n",
	  TSTART, TEND, DT);

  FILE* monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fprintf (monitor, "\n#####################\n");
  fprintf (monitor, "Program fFileGenerate\n");
  fprintf (monitor, "#####################\n");
  fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  
  fprintf (monitor, "Input parameters (from Inputs/fFile.nml):\n");
  fprintf (monitor, "Git Hash     = "); fprintf (monitor, GIT_HASH);     fprintf (monitor, "\n");
  fprintf (monitor, "Compile time = "); fprintf (monitor, COMPILE_TIME); fprintf (monitor, "\n");
  fprintf (monitor, "Git Branch   = "); fprintf (monitor, GIT_BRANCH);   fprintf (monitor, "\n\n");
  fprintf (monitor, "FLUX_NTOR = %2d  FLUX_MMIN  = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  fprintf (monitor, "NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	   NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  fprintf (monitor, "TSTART = %11.4e  TEND = %11.4e  DT = %11.4e\n",
	  TSTART, TEND, DT);
  fclose (monitor);

  FILE* namelist = fopen ("Inputs/InputParameters", "w");
  fprintf (namelist, "Input parameters (from Inputs/fFile.nml):\n");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "FLUX_NTOR = %2d  FLUX_MMIN = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  fprintf (namelist, "NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	   NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  fprintf (namelist, "TSTART = %11.4e  TEND = %11.4e  DT = %11.4e\n",
	  TSTART, TEND, DT);
  fclose (namelist);

  // ............
  // Perform scan
  // ............
  double Time = TSTART;
  char   FLUXstring[MAXCOMMANDLINELENGTH], NEOstring[MAXCOMMANDLINELENGTH];

  do
    {
      printf  ("\n$$$$$$$$$$$$$$$$$$\ntime = %11.4e\n$$$$$$$$$$$$$$$$$$\n", Time);

      monitor = fopen ("Outputs/monitor.txt", "a");
      fprintf  (monitor, "\n$$$$$$$$$$$$$$$$$$\ntime = %11.4e\n$$$$$$$$$$$$$$$$$$\n", Time);
      fclose (monitor);

      // Construct command strings
      sprintf (FLUXstring,  "cd ../Flux; ./flux -g %d -n %d -m %d -M %d -t %16.9e ",
	       FLUX_INTG, FLUX_NTOR, FLUX_MMIN, FLUX_MMAX, Time);
      sprintf (NEOstring,   "cd ../Neoclassical; ./neoclassical -e %d -p %d -n %d -I %d -f %2d -T %2d -N %16.9e -l %16.9e -y %16.9e -t %16.9e ",
	       NEO_INTF, NEO_INTP, NEO_NEUTRAL, NEO_IMPURITY, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN, Time);
  
      // Call program FLUX
      printf ("Executing:: %s\n", FLUXstring);
      if (CallSystem (FLUXstring) != 0)
	exit (1);
      
      // Call program NEOCLASSICAL
      printf ("Executing:: %s\n", NEOstring);
      if (CallSystem (NEOstring) != 0)
	exit (1);

      // Increment time
      Time += DT;
    }
  while (Time <= TEND);

  time                 (&rawtime);
  timeinfo = localtime (&rawtime);

  monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fclose (monitor);
}