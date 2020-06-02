// main.cpp

// ###########################################################
// Program to simulate mulit-harmomic magnetic island dynamics
// in presence of RMP in time-evolving toroidal tokamak
// discharge.
//
// Uses FLUX/NEOCLASSICAL/PHASE codes.
// ###########################################################

#include <stdio.h>
#include <stdlib.h>

extern "C" void NameListRead (int* IMPURITY, int* NEUTRAL, int* FREQ, int* STAGE2,
			      int* INTP, int* RESTART, int* OLD, double* TSTART, double* TEND, double* DT);

// #############
// Main function
// #############

int main ()
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n#####################\n");
  printf ("Program IslandDyamics\n");
  printf ("######################\n\n");

  // .............
  // Read namelist
  // .............
  int    IMPURITY; // If != 0 then single ion impurity species included in calculation
  int    NEUTRAL;  // If != 0 then majority ion neutrals included in calculation
  int    FREQ;     // If < 0 || == 0 || > 0 then use linear/nonlinear/ExB natural frequency
  int    INTP;     // If == 0 then fFile/nfile/uFile/lFile used
	           // If == 1 then fFile/nFile interpolated and uFile/lFile used
		   // If == 2 then all files interpolated
  int    STAGE2;   // If != 0 then Stage2 PHASE calculation enabled
  int    RESTART;  // If != 0 then delete all previous EPEC data
  int    OLD;      // If != 0 then restart PHASE calculations from previous run
  double TSTART;   // Simulation experimental start time (ms)
  double TEND;     // Simulation experimental end time (ms)
  double DT;       // Simulation experimental time step (ms)

  NameListRead (&IMPURITY, &NEUTRAL, &FREQ, &STAGE2, &INTP, &RESTART, &OLD, &TSTART, &TEND, &DT);

  printf ("Reading namelist:\n");
  printf ("IMPURITY = %2d  NEUTRAL = %2d  FREQ = %2d  INTP = %2d  STAGE2 = %2d  RESTART = %2d  OLD = %2d  TSTART = %11.4e  TEND = %11.4e  DT = %11.4e\n",
	  IMPURITY, NEUTRAL, FREQ, INTP, STAGE2, RESTART, OLD, TSTART, TEND, DT);
 
  // Sanity check
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

  // Initialize symbolic links
  system ("cd ../Flux/fFiles; rm -rf *");
  system ("cd ../Neoclassical/nFiles; rm -rf *");

  // Initialize output files
  if (RESTART)
    system ("cd Stage1; rm *.txt");
  
  // ..................
  // Perform simulation
  // ..................
  double time = TSTART;
  char   FLUXstring[200], NEOstring[200], PHASEstring[200];
  
  FILE* file = fopen ("Stage1/monitor.txt", "w");
  
  do
    {
      printf  ("\ntime = %11.4e\n", time);

      // Update monitor file
      fprintf (file, "time = %11.4e\n", time);
      fflush  (file);

      // Construct command strings
      sprintf (FLUXstring,  "cd ../Flux; ./flux -i 1 -t %16.9e ",
	       time);
      sprintf (NEOstring,   "cd ../Neoclassical; ./neoclassical -i 1 -n %2d -I %2d -f %2d -t %16.9e ",
	       NEUTRAL, IMPURITY, FREQ, time);
      sprintf (PHASEstring, "cd ../Phase; ./phase -i %d -s %2d -o %2d -t %16.9e ",
	       INTP, STAGE2, OLD, time);
      
      // Call program FLUX
      if (INTP < 2)
	if (system (FLUXstring) != 0)
	  exit (1);

      // Call program NEOCLASSICAL
      if (INTP < 2)
	if (system (NEOstring) != 0)
	  exit (1);

      // Call program PHASE
      if (system (PHASEstring) != 0)
	exit (1);

      // Increment time
      time += DT;
    }
  while (time < TEND);

  fclose (file);
}
