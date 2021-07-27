// IslandDynamics.cpp

// ##############################################################
// Program to simulate multi-harmonic magnetic island dynamics
// in presence of static, externally generated, resonant magnetic
// perturbation in time-evolving toroidal tokamak discharge.

// Uses FLUX/NEOCLASSICAL/PHASE codes.

// Namelist in Inputs/Island.nml

// .........
// Versions:
// .........

// 1.0  - Initial version
// 1.1  - Added time to monitor.txt
// 1.2  - Major rearrangement of input and output files
// 1.3  - Added PHASE_FREQ flag
// 1.4  - Added PHASE_SCALE flag
// 1.5  - Added PHASE_LIN flag
// 1.6  - Added PHASE_MID flag
// 1.7  - Added PHASE_CHIR flag
// 1.8  - Added IRMP
// 1.9  - Added COPT flag
// 1.10 - Added more copt options
// 1.11 - Added PSIPED
// 1.12 - Removed PHASE_FREQ flag
// 1.13 - Greatly simplified such that only calls PHASE
// 1.14 - Added PHASE_CXD and PHASE_BSC flags
// 1.15 - Added PHASE_BOOT, PHASE_CURV, PHASE_POLZ flags

// ##############################################################

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 14

#define MAXCOMMANDLINELENGTH 500

extern "C" void NameListRead (int* PHASE_MID, int* PHASE_COPT, double* PHASE_CORE, int* PHASE_LIN, int* PHASE_FREQ,
			      double* PHASE_FFAC, int* PHASE_HIGH, int* PHASE_RATS, int* PHASE_NATS, double* PHASE_SCALE,
			      double* PHASE_CHIR, int* PHASE_CXD, int* PHASE_BOOT, int* PHASE_CURV, int* PHASE_POLZ,
			      int* RESTART, double* TSTART, double* TEND, double* DT); 

void IslandDynamics (double IRMP);

// #############
// Main function
// #############

int main (int argc, char** argv)
{
  int c;
  char* ivalue = NULL;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "hi:")) != -1)
    switch (c)
      {
      case 'h':
 	printf ("Options:\n");
 	printf ("-h      - list options\n");
	printf ("-i IRMP - set RMP current to IRMP");
	exit (0);
      case 'i':
	ivalue = optarg;
 	break;
      case '?':
	if (optopt == 'i')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  double _IRMP = -1.e6; float __IRMP;
  
  if (ivalue != NULL)
     {
      __IRMP = atof (ivalue);
      _IRMP  = double (__IRMP);
     }
    
  IslandDynamics (_IRMP);

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

// ##############################################
// Function to perform island dynamics simulation
// ##############################################
void IslandDynamics (double IRMP)
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n######################\n");
  printf ("Program IslandDynamics\n");
  printf ("######################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // .............
  // Read namelist
  // .............
  int	 PHASE_STAGE5 = 1;   // If != 0 then Stage5 PHASE calculation enabled
  int	 PHASE_INTF   = 1;   // If != 0 then interpolate fFiles in PHASE
  int	 PHASE_INTN   = 1;   // If != 0 then interpolate nFiles in PHASE
  int	 PHASE_INTU   = 1;   // If != 0 then interpolate uFiles/mFiles/lFiles in PHASE
  int	 PHASE_OLD    = 1;   // If != 0 then restart PHASE calculations from previous run

  int	 PHASE_MID;          // If != 0 then middle RMP coils included in PHASE calculation
  int    PHASE_COPT;         // If != 0 then optimize coil currents
  double PHASE_CORE;         // Core drive minimization factor

  int    PHASE_LIN;          // If != 0 then perform purely linear calculation in PHASE
  int    PHASE_FREQ;         // Flag for selecting natural frequency type in PHASE
                             //  If == 0 then use composite linear/nonlinear model
                             //  If == 1 then w_natural = FFAC * w_linear + (1-FFAC) * w_EB
  double PHASE_FFAC;         // Parameter for determining natural frequency in PHASE

  int    PHASE_CXD;          // If != 0 then include charge exchange damping in angular equations of motion in PHASE
  int    PHASE_BOOT;         // If != 0 then include bootstrap terms in Rutherford equations in PHASE
  int    PHASE_CURV;         // If != 0 then include curvature terms in Rutherford equations in PHASE
  int    PHASE_POLZ;         // If != 0 then include polarization terms in Rutherford equations in PHASE

  double PHASE_SCALE;        // GPEC scalefactor in PHASE
  double PHASE_CHIR;         // Maximum Chirikov parameter for vacuum islands in PHASE
  int    PHASE_HIGH;         // If != 0 then perform higher order transport calculation in PHASE
  int    PHASE_RATS;         // If != 0 then perform linear-only interpolation of uFiles/mFiles/lFiles in PHASE
  int    PHASE_NATS;         // If != 0 then perform linear-only interpolation of nFiles in PHASE

  int	 RESTART;            // If != 0 then delete all previous IslandDynamics data
  double TSTART;             // Simulation experimental start time (ms)
  double TEND;               // Simulation experimental end time (ms)
  double DT;                 // Simulation experimental time step (ms)

  NameListRead (&PHASE_MID, &PHASE_COPT, &PHASE_CORE, &PHASE_LIN, &PHASE_FREQ, &PHASE_FFAC, 
		&PHASE_HIGH, &PHASE_RATS, &PHASE_NATS, &PHASE_SCALE, &PHASE_CHIR, &PHASE_CXD,
		&PHASE_BOOT, &PHASE_CURV, &PHASE_POLZ,
		&RESTART, &TSTART, &TEND, &DT);

  printf ("Reading Inputs/Island.nml:\n");
  printf ("PHASE_MID = %2d PHASE_COPT = %2d PHASE_CORE = %11.4e PHASE_LIN = %2d PHASE_FREQ = %2d PHASE_FFAC = %11.4e PHASE_HIGH = %2d PHASE_RATS = %2d PHASE_NATS = %2d PHASE_SCALE = %11.4e PHASE_CHIR = %11.4e PHASE_CXD = %2d PHASE_BSC = %2d PHASE_CURV = %2d PHASE_POLZ = %2d\d\n",
	  PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_LIN, PHASE_FREQ, PHASE_FFAC, PHASE_HIGH, PHASE_RATS, PHASE_NATS, PHASE_SCALE, PHASE_CHIR, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ);
  printf ("RESTART = %2d TSTART = %11.4e TEND = %11.4e DT = %11.4e IRMP = %11.4e\n",
	  RESTART, TSTART, TEND, DT, IRMP);

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

  // .......................
  // Initialize output files
  // .......................
  if (RESTART)
    {
      CallSystem ("cd Outputs/Stage6; rm *.txt");
      CallSystem ("cd Outputs; rm monitor.txt");
    }

  // .................
  // Get date and time
  // .................
  time_t     rawtime;
  struct tm* timeinfo;
  time                 (&rawtime);
  timeinfo = localtime (&rawtime);
  
  printf ("%s\n", asctime (timeinfo));
  printf ("\n######################\n");
  printf ("Program IslandDynamics\n");
  printf ("######################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  
  printf ("Input parameters (from Inputs/Island.nml):\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("PHASE_MID = %2d PHASE_COPT = %2d PHASE_CORE = %11.4e PHASE_LIN = %2d PHASE_FREQ = %2d PHASE_FFAC = %11.4e PHASE_HIGH = %2d PHASE_RATS = %2d PHASE_NATS = %2d PHASE_SCALE = %11.4e PHASE_CHIR = %11.4e PHASE_CXD = %2d PHASE_BOOT = %2d PHASE_CURV = %2d PHASE_POLZ = %2d\n",
	  PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_LIN, PHASE_FREQ, PHASE_FFAC, PHASE_HIGH, PHASE_RATS, PHASE_NATS, PHASE_SCALE, PHASE_CHIR, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ);
  printf ("RESTART = %2d TSTART = %11.4e TEND = %11.4e DT = %11.4e IRMP = %11.4e\n",
	  RESTART, TSTART, TEND, DT, IRMP);

   FILE* monitor = fopen ("Outputs/monitor.txt", "a");
  fprintf (monitor, "%s\n", asctime (timeinfo));
  fprintf (monitor, "\n######################\n");
  fprintf (monitor, "Program IslandDynamics\n");
  fprintf (monitor, "######################\n");
  fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  
  fprintf (monitor, "Input parameters (from Inputs/Island.nml):\n");
  fprintf (monitor, "Git Hash     = "); fprintf (monitor, GIT_HASH);     fprintf (monitor, "\n");
  fprintf (monitor, "Compile time = "); fprintf (monitor, COMPILE_TIME); fprintf (monitor, "\n");
  fprintf (monitor, "Git Branch   = "); fprintf (monitor, GIT_BRANCH);   fprintf (monitor, "\n\n");
  fprintf (monitor, "PHASE_MID = %2d PHASE_COPT = %2d PHASE_CORE = %11.4e PHASE_LIN = %2d PHASE_FREQ = %2d PHASE_FFAC = %11.4e PHASE_HIGH = %2d PHASE_RATS = %2d PHASE_NATS = %2d PHASE_SCALE = %11.4e PHASE_CHIR = %11.4e PHASE_CXD = %2d PHASE_BOOT = %2d PHASE_CURV = %2d PHASE_POLZ = %2d\n",
	   PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_LIN, PHASE_FREQ, PHASE_FFAC, PHASE_HIGH, PHASE_RATS, PHASE_NATS, PHASE_SCALE, PHASE_CHIR, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ);
  fprintf (monitor, "RESTART = %2d TSTART = %11.4e TEND = %11.4e DT = %11.4e IRMP = %11.4e\n",
	  RESTART, TSTART, TEND, DT, IRMP);
  fclose (monitor);

  FILE* namelist = fopen ("Inputs/InputParameters", "w");
  fprintf (namelist, "Input parameters (from Inputs/Island.nml):\n");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "PHASE_MID = %2d PHASE_COPT = %2d PHASE_CORE = %11.4e PHASE_LIN = %2d PHASE_FREQ = %2d PHASE_FFAC = %11.4e PHASE_HIGH = %2d PHASE_RATS = %2d PHASE_NATS = %2d PHASE_SCALE = %11.4e PHASE_CHIR = %11.4e PHASE_CXD = %2d PHASE_BOOT = %2d PHASE_CURV = %2d PHASE_POLZ = %2d\n",
	   PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_LIN, PHASE_FREQ, PHASE_FFAC, PHASE_HIGH, PHASE_RATS, PHASE_NATS, PHASE_SCALE, PHASE_CHIR, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ);
  fprintf (namelist, "RESTART = %2d TSTART = %11.4e TEND = %11.4e DT = %11.4e IRMP = %11.4e\n",
	  RESTART, TSTART, TEND, DT, IRMP);
  fclose (namelist);

  // ..................
  // Perform simulation
  // ..................
  double Time = TSTART;
  char   PHASEstring[MAXCOMMANDLINELENGTH];

  int lock = 0;
  do
    {
      printf  ("\n$$$$$$$$$$$$$$$$$$\ntime = %11.4e\n$$$$$$$$$$$$$$$$$$\n", Time);

      monitor = fopen ("Outputs/monitor.txt", "a");
      fprintf  (monitor, "\n$$$$$$$$$$$$$$$$$$\ntime = %11.4e\n$$$$$$$$$$$$$$$$$$\n", Time);
      fclose (monitor);
  
      if (lock)
	sprintf (PHASEstring, "cd ../Phase; ./phase -m %d -C %d -f %d -n %d -u %d -s %2d -o %2d -l %2d -t %16.9e -T %16.9e -S %16.9e -c %11.4e -i %11.4e -H %2d -r %2d -D %16.9e -F %2d -N %2d -G %16.9e -X %2d -B %2d -V %2d -P %2d",
		 PHASE_MID, PHASE_COPT, PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_LIN, Time, Time + DT,
		 PHASE_SCALE, PHASE_CHIR, IRMP, PHASE_HIGH, PHASE_RATS, PHASE_CORE, PHASE_FREQ, PHASE_NATS, PHASE_FFAC, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ);
      else
	sprintf (PHASEstring, "cd ../Phase; ./phase -m %d -C %d -f %d -n %d -u %d -s %2d -o %2d -l %2d -t %16.9e -T %16.9e -S %16.9e -c %11.4e -i %11.4e -H %2d -r %2d -D %16.9e -F %2d -N %2d -G %16.9e -X %2d -B %2d -V %2d -P %2d",
		 PHASE_MID, PHASE_COPT, PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, lock, PHASE_LIN, Time, Time + DT,
		 PHASE_SCALE, PHASE_CHIR, IRMP, PHASE_HIGH, PHASE_RATS, PHASE_CORE, PHASE_FREQ, PHASE_NATS, PHASE_FFAC, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ);
      
      lock = 1;
      
      // Call PHASE
      printf ("Executing:: %s\n", PHASEstring);
      if (CallSystem (PHASEstring) != 0)
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
