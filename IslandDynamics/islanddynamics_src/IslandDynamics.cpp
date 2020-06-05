// IslandDynamics.cpp

// ##############################################################
// Program to simulate multi-harmonic magnetic island dynamics
// in presence of static, externally generated, resonant magnetic
// perturbation in time-evolving toroidal tokamak discharge.
//
// Uses FLUX/NEOCLASSICAL/PHASE codes.
// ##############################################################

#include <stdio.h>
#include <stdlib.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

#define MAXCOMMANDLINELENGTH 500

extern "C" void NameListRead (int* FLUX_NTOR, int* FLUX_MMIN, int* FLUX_MMAX, int* NEO_INTF,
			      int* NEO_IMPURITY, int* NEO_NEUTRAL, int* NEO_FREQ, double* NEO_YN,
			      int* PHASE_INTN, int* PHASE_STAGE2, int* PHASE_OLD,
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
  printf ("#####################\n\n");
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
  double NEO_YN;         // Set neutral poloidal asymmetry parameter in NEOCLASSICAL
  int	 PHASE_INTF;     // If != 0 then interpolate fFiles in PHASE
  int	 PHASE_INTN;     // If != 0 then interpolate nFiles in PHASE
  int	 PHASE_INTU = 1; // If != 0 then interpolate uFiles/lFiles in PHASE
  int	 PHASE_STAGE2;   // If != 0 then Stage2 PHASE calculation enabled
  int	 PHASE_OLD;      // If != 0 then restart PHASE calculations from previous run
  int	 RESTART;        // If != 0 then delete all previous IslandDynamics data
  double TSTART;         // Simulation experimental start time (ms)
  double TEND;           // Simulation experimental end time (ms)
  double DT;             // Simulation experimental time step (ms)

  NameListRead (&FLUX_NTOR, &FLUX_MMIN, &FLUX_MMAX, &NEO_INTF, &NEO_IMPURITY, &NEO_NEUTRAL, &NEO_FREQ, &NEO_YN,
		&PHASE_INTN, &PHASE_STAGE2, &PHASE_OLD, &RESTART, &TSTART, &TEND, &DT);

  PHASE_INTF = NEO_INTF;

  printf ("Reading namelist:\n");
  printf ("FLUX_NTOR  = %2d  FLUX_MMIN    = %2d           FLUX_MMAX   = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  printf ("NEO_INTF   = %2d  NEO_IMPURITY = %2d           NEO_NEUTRAL = %2d           NEO_FREQ = %2d  NEO_YN = %11.4e\n",
	  NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_YN);
  printf ("PHASE_INTN = %2d  PHASE_STAGE2 = %2d           PHASE_OLD   = %2d\n",
	  PHASE_INTN, PHASE_STAGE2, PHASE_OLD);
  printf ("RESTART    = %2d  TSTART       = %11.4e  TEND        = %11.4e  DT       = %11.4e\n",
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
  system ("cd ../Flux/fFiles; rm -rf *");
  system ("cd ../Neoclassical/nFiles; rm -rf *");

  // .......................
  // Initialize output files
  // .......................
  if (RESTART)
    system ("cd Stage1; rm *.txt");
  
  // ..................
  // Perform simulation
  // ..................
  double time = TSTART;
  char   FLUXstring[MAXCOMMANDLINELENGTH], NEOstring[MAXCOMMANDLINELENGTH], PHASEstring[MAXCOMMANDLINELENGTH];
  
  FILE* file = fopen ("Stage1/monitor.txt", "w");
  
  do
    {
      printf  ("\n$$$$$$$$$$$$$$$$$$\n$$$$$$$$$$$$$$$$$$\ntime = %11.4e\n$$$$$$$$$$$$$$$$$$\n$$$$$$$$$$$$$$$$$$\n\n", time);

      // Update monitor file
      fprintf (file, "time = %11.4e\n", time);
      fflush  (file);

      // Construct command strings
      sprintf (FLUXstring,  "cd ../Flux; ./flux -g %d -n %d -m %d -M %d -t %16.9e ",
	       FLUX_INTG, FLUX_NTOR, FLUX_MMIN, FLUX_MMAX, time);
      sprintf (NEOstring,   "cd ../Neoclassical; ./neoclassical -e %d -i %d -n %d -I %d -f %2d -y %16.9e -t %16.9e ",
	       NEO_INTF, NEO_INTP, NEO_NEUTRAL, NEO_IMPURITY, NEO_FREQ, NEO_YN, time);
      sprintf (PHASEstring, "cd ../Phase; ./phase -f %d -n %d -u %d -s %2d -o %2d -t %16.9e ",
	       PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE2, PHASE_OLD, time);
      
      // Call program FLUX
      if (NEO_INTF == 0)
	{
	  printf ("%s\n", FLUXstring);
	  if (system (FLUXstring) != 0)
	    exit (1);
	}

      // Call program NEOCLASSICAL
      if (PHASE_INTN == 0)
	{
	  printf ("%s\n", NEOstring);
	  if (system (NEOstring) != 0)
	    exit (1);
	}

      // Call program PHASE
      printf ("%s\n", PHASEstring);
      if (system (PHASEstring) != 0)
	exit (1);

      // Increment time
      time += DT;
    }
  while (time < TEND);

  fclose (file);
}
