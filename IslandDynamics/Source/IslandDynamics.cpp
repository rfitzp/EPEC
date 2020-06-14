// IslandDynamics.cpp

// ##############################################################
// Program to simulate multi-harmonic magnetic island dynamics
// in presence of static, externally generated, resonant magnetic
// perturbation in time-evolving toroidal tokamak discharge.

// Uses FLUX/NEOCLASSICAL/PHASE codes.

// Version:

// 1.0 - Initial version
// 1.1 - Added time to monitor.txt
// 1.2 - Major rearrangement of input and output files

// ##############################################################

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 2

#define MAXCOMMANDLINELENGTH 500

extern "C" void NameListRead (int* FLUX_NTOR, int* FLUX_MMIN, int* FLUX_MMAX,
			      int* NEO_INTF, int* NEO_IMPURITY, int* NEO_NEUTRAL, int* NEO_FREQ, int* NEO_NTYPE, double* NEO_NN, double* NEO_LN, double* NEO_YN,
			      int* PHASE, int* PHASE_INTN, int* PHASE_STAGE2, int* PHASE_OLD,
			      int* RESTART, double* TSTART, double* TEND, double* DT); 

void IslandDynamics ();

// #############
// Main function
// #############

int main ()
{
  IslandDynamics ();

  return 0;
}

// ##############################################
// Function to perform island dynamics simulation
// ##############################################
void IslandDynamics ()
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n#####################\n");
  printf ("Program IslandDyamics\n");
  printf ("#####################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // .............
  // Read namelist
  // .............
  int    FLUX_INTG = 1;  // If != 0 then interpolate gFiles in FLUX
  int    FLUX_NTOR;      // Set toroidal mode number in FLUX
  int    FLUX_MMIN;      // Set minimum poloidal mode number in FLUX
  int	 FLUX_MMAX;      // Set maximum poloidal mode number in FLUX
  int    NEO_INTP = 1;   // If != 0 then interpolate pFiles in NEOCLASSICAL
  int	 NEO_INTF;       // If != 0 then interpolate fFiles in NEOCLASSICAL
  int	 NEO_IMPURITY;   // If != 0 then include single ion impurity species in NEOCLASSICAL
  int	 NEO_NEUTRAL;    // If != 0  then include majority ion neutrals in NEOCLASSICAL
  int	 NEO_FREQ;       // If < 0 || == 0 || > 0 then use linear/nonlinear/ExB natural frequency in NEOCLASSICAL
  int    NEO_NTYPE;      // If 0/1 then neutral density distribution exponential/Lorentzian in NEOCLASSICAL
  double NEO_NN;         // Set neutral density at boundary in NEOCLASSICAL
  double NEO_LN;         // Set neutral density scalelength in NEOCLASSICAL
  double NEO_YN;         // Set neutral poloidal asymmetry parameter in NEOCLASSICAL
  int    PHASE;          // If != 0 then call PHASE
  int	 PHASE_INTF;     // If != 0 then interpolate fFiles in PHASE
  int	 PHASE_INTN;     // If != 0 then interpolate nFiles in PHASE
  int	 PHASE_INTU = 1; // If != 0 then interpolate uFiles/lFiles in PHASE
  int	 PHASE_STAGE5;   // If != 0 then Stage5 PHASE calculation enabled
  int	 PHASE_OLD;      // If != 0 then restart PHASE calculations from previous run
  int	 RESTART;        // If != 0 then delete all previous IslandDynamics data
  double TSTART;         // Simulation experimental start time (ms)
  double TEND;           // Simulation experimental end time (ms)
  double DT;             // Simulation experimental time step (ms)

  NameListRead (&FLUX_NTOR, &FLUX_MMIN, &FLUX_MMAX,
		&NEO_INTF, &NEO_IMPURITY, &NEO_NEUTRAL, &NEO_FREQ, &NEO_NTYPE, &NEO_NN, &NEO_LN, &NEO_YN,
		&PHASE, &PHASE_INTN, &PHASE_STAGE5, &PHASE_OLD,
		&RESTART, &TSTART, &TEND, &DT);

  PHASE_INTF = NEO_INTF;

  printf ("Reading Inputs/Island.in:\n");
  printf ("FLUX_NTOR  = %2d  FLUX_MMIN    = %2d           FLUX_MMAX    = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  printf ("NEO_INTF   = %2d  NEO_IMPURITY = %2d           NEO_NEUTRAL  = %2d           NEO_FREQ  = %2d           NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  printf ("PHASE      = %2d  PHASE_INTN   = %2d           PHASE_STAGE5 = %2d           PHASE_OLD = %2d\n",
	  PHASE, PHASE_INTN, PHASE_STAGE5, PHASE_OLD);
  printf ("RESTART    = %2d  TSTART       = %11.4e  TEND         = %11.4e  DT        = %11.4e\n",
	  RESTART, TSTART, TEND, DT);

  // ............
  // Sanity check
  // ............
  if (TEND < TSTART)
    {
      printf ("Error - TEND cannot be less than TSTART\n");
      exit (1);
    }
  if (DT <= 0.)
    {
      printf ("Error - DT must be positive\n");
      exit (1);
    }

  // .........................
  // Initialize symbolic links
  // .........................
  system ("cd ../Flux/Outputs/fFiles; rm -rf *");
  system ("cd ../Neoclassical/Outputs/nFiles; rm -rf *");

  // .......................
  // Initialize output files
  // .......................
  if (RESTART)
    {
      system ("cd Outputs/Stage6; rm *.txt");
      system ("cd Outputs; rm monitor.txt");
    }

  // .................
  // Get date and time
  // .................
  time_t      rawtime;
  struct tm * timeinfo;
  time                 (&rawtime);
  timeinfo = localtime (&rawtime);
 
  FILE* monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fprintf (monitor, "\n######################\n");
  fprintf (monitor, "Program IslandDynamics\n");
  fprintf (monitor, "######################\n");
  fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  fprintf (monitor, "Input parameters (from Inputs/Island.in):\n");
  fprintf (monitor, "FLUX_NTOR  = %2d  FLUX_MMIN    = %2d           FLUX_MMAX    = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  fprintf (monitor, "NEO_INTF   = %2d  NEO_IMPURITY = %2d           NEO_NEUTRAL  = %2d           NEO_FREQ  = %2d           NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  fprintf (monitor, "PHASE      = %2d  PHASE_INTN   = %2d           PHASE_STAGE5 = %2d           PHASE_OLD = %2d\n",
	  PHASE, PHASE_INTN, PHASE_STAGE5, PHASE_OLD);
  fprintf (monitor, "RESTART    = %2d  TSTART       = %11.4e  TEND         = %11.4e  DT        = %11.4e\n",
	  RESTART, TSTART, TEND, DT);
  fclose (monitor);

  FILE* namelist = fopen ("Inputs/InputParameters", "w");
  fprintf (namelist, "Input parameters (from Inputs/Island.in):\n");
  fprintf (namelist, "FLUX_NTOR  = %2d  FLUX_MMIN    = %2d           FLUX_MMAX    = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  fprintf (namelist, "NEO_INTF   = %2d  NEO_IMPURITY = %2d           NEO_NEUTRAL  = %2d           NEO_FREQ  = %2d           NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  fprintf (namelist, "PHASE      = %2d  PHASE_INTN   = %2d           PHASE_STAGE5 = %2d           PHASE_OLD = %2d\n",
	  PHASE, PHASE_INTN, PHASE_STAGE5, PHASE_OLD);
  fprintf (namelist, "RESTART    = %2d  TSTART       = %11.4e  TEND         = %11.4e  DT        = %11.4e\n",
	  RESTART, TSTART, TEND, DT);
  fclose (namelist);

  // ..................
  // Perform simulation
  // ..................
  double Time = TSTART;
  char   FLUXstring[MAXCOMMANDLINELENGTH], NEOstring[MAXCOMMANDLINELENGTH], PHASEstring[MAXCOMMANDLINELENGTH];

  int lock = 0;
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
      if (lock)
	sprintf (PHASEstring, "cd ../Phase; ./phase -f %d -n %d -u %d -s %2d -o %2d -t %16.9e ",
		 PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, Time);
      else
	sprintf (PHASEstring, "cd ../Phase; ./phase -f %d -n %d -u %d -s %2d -o %2d -t %16.9e ",
		 PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, lock, Time);
      lock = 1;
      
      // Call program FLUX
      if (NEO_INTF == 0)
	{
	  printf ("Executing:: %s\n", FLUXstring);
	  if (system (FLUXstring) != 0)
	    exit (1);
	}

      // Call program NEOCLASSICAL
      if (PHASE_INTN == 0)
	{
	  printf ("Executing:: %s\n", NEOstring);
	  if (system (NEOstring) != 0)
	    exit (1);
	}

      // Call program PHASE
      if (PHASE)
	{
	  printf ("Executing:: %s\n", PHASEstring);
	  if (system (PHASEstring) != 0)
	    exit (1);
	}

      // Increment time
      Time += DT;
    }
  while (Time < TEND);

  time                 (&rawtime);
  timeinfo = localtime (&rawtime);

  monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fclose (monitor);
}
