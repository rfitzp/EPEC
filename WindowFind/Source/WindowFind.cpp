// WindowFind.cpp

// ##############################################################
// Program to find RMP-induced ELM suppression windows in q95
// versus Irmp space.
// Uses FLUX/NEOCLASSICAL/PHASE/ISLANDDYNAMICS codes.

// .........
// Versions:
// .........

// 1.0 - Initial version

// ##############################################################

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

#define MAXCOMMANDLINELENGTH 500

extern "C" void NameListRead (double* ISTART, double* IEND, double* DI); 

void WindowFind ();

// #############
// Main function
// #############

int main (int argc, char** argv)
{
  WindowFind ();

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

// ########################################
// Function to find EKM suppression windows
// ########################################
void WindowFind ()
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n##################\n");
  printf ("Program WindowFind\n");
  printf ("##################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // .............
  // Read namelist
  // .............
  double ISTART;  // Initial RMP current time (kA)
  double IEND;    // Final RMP current (kA)
  double DI;      // RMP current step (kA)

  NameListRead (&ISTART, &IEND, &DI);

  printf ("Reading Inputs/Window.nml:\n");
  printf ("ISTART = %11.4e  IEND = %11.4e  DI = %11.4e\n",
	  ISTART, IEND, DI);

  // ............
  // Sanity check
  // ............
  if (IEND < ISTART)
    {
      printf ("Error - IEND must be greater than ISTART\n");
      exit (1);
    }
  if (DI < 0.)
    {
      printf ("Error - DI must be positive\n");
      exit (1);
    }

  // .........................
  // Initialize symbolic links
  // .........................
  CallSystem ("cd ../Flux/Outputs/fFiles; rm -rf *");
  CallSystem ("cd ../Neoclassical/Outputs/nFiles; rm -rf *");

  // .......................
  // Initialize output files
  // .......................
  CallSystem ("cd Outputs; rm *.txt");
  CallSystem ("cd Outputs; rm monitor.txt");
  
  // .................
  // Get date and time
  // .................
  time_t     rawtime;
  struct tm* timeinfo;
  time                 (&rawtime);
  timeinfo = localtime (&rawtime);
  
  printf ("%s\n", asctime (timeinfo));
  printf ("\n##################\n");
  printf ("Program WindowFind\n");
  printf ("##################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  
  printf ("Input parameters (from Inputs/Window.nml):\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("ISTART = %11.4e  IEND = %11.4e  DI = %11.4e\n",
	  ISTART, IEND, DI);

  FILE* monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fprintf (monitor, "\n##################\n");
  fprintf (monitor, "Program WindowFind\n");
  fprintf (monitor, "##################\n");
  fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  
  fprintf (monitor, "Input parameters (from Inputs/Window.nml):\n");
  fprintf (monitor, "Git Hash     = "); fprintf (monitor, GIT_HASH);     fprintf (monitor, "\n");
  fprintf (monitor, "Compile time = "); fprintf (monitor, COMPILE_TIME); fprintf (monitor, "\n");
  fprintf (monitor, "Git Branch   = "); fprintf (monitor, GIT_BRANCH);   fprintf (monitor, "\n\n");
  fprintf (monitor, "ISTART = %11.4e  IEND = %11.4e  DI = %11.4e\n",
	   ISTART, IEND, DI);
  fclose (monitor);

  FILE* namelist = fopen ("Inputs/InputParameters.txt", "w");
  fprintf (namelist, "Input parameters (from Inputs/Window.nml):\n");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "ISTART = %11.4e  IEND = %11.4e  DI = %11.4e\n",
	   ISTART, IEND, DI);
  fclose (namelist);

  // ..............
  // Perform search
  // ..............
  double IRMP = ISTART;
  char   ISLANDstring[MAXCOMMANDLINELENGTH];
  double a, b, c, Q[5000], P[5000];

  do
    {
      printf  ("\n$$$$$$$$$$$$$$$$$$\nIRMP = %11.4e\n$$$$$$$$$$$$$$$$$$\n", IRMP);

      monitor = fopen ("Outputs/monitor.txt", "a");
      fprintf  (monitor, "\n$$$$$$$$$$$$$$$$$$\nIRMP = %11.4e\n$$$$$$$$$$$$$$$$$$\n", IRMP);
      fclose (monitor);

      // Construct command strings
      sprintf (ISLANDstring,  "cd ../IslandDynamics; ./island -i %e ", IRMP);
          
      // Call program ISLANDDYNAMICS
      printf ("Executing:: %s\n", ISLANDstring);
      if (CallSystem (ISLANDstring) != 0)
	exit (1);

      // Extract RMP window data
      FILE* file = fopen ("../IslandDynamics/Outputs/Stage6/deltap.txt", "r");
      int i = -1;
      while (fscanf (file, "%lf %lf %lf", &a, &b, &c) == 3)
	{
	  i++;
	  Q[i] = b;
	  P[i] = c;
	}
      fclose (file);

      // Write RMP window data
      if (Q[0] < Q[i])
	{
	  file = fopen ("Outputs/window.txt", "a");
  	  for (int j = 0; j <= i; j++)
	    fprintf (file, "%16.9e ", P[j]);
	  fprintf (file, "\n");
	  fclose (file);
	  
	  file = fopen ("Outputs/limits.txt", "w");
	  fprintf (file, "%16.9e %16.9e %16.9e %16.9e\n", Q[0], Q[i], ISTART, IRMP);
	  fclose (file);
	}
      else
	{
	  file = fopen ("Outputs/window.txt", "a");
  	  for (int j = 0; j <= i; j++)
	    fprintf (file, "%16.9e ", P[i-j]);
	  fprintf (file, "\n");
	  fclose (file);
	  
	  file = fopen ("Outputs/limits.txt", "w");
	  fprintf (file, "%16.9e %16.9e %16.9e %16.9e\n", Q[i], Q[0], ISTART, IRMP);
	  fclose (file);
	}

      // Increment RMP current
      IRMP += DI;
    }
  while (IRMP <= IEND);

  time                 (&rawtime);
  timeinfo = localtime (&rawtime);

  monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fclose (monitor);
}

