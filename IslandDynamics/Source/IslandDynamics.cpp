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

// ##############################################################

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 10

#define MAXCOMMANDLINELENGTH 500

extern "C" void NameListRead (int* FLUX_NTOR, int* FLUX_MMIN, int* FLUX_MMAX,
			      int* NEO_INTP, int* NEO_INTF, int* NEO_IMPURITY, int* NEO_NEUTRAL, int* NEO_FREQ, int* NEO_NTYPE, double* NEO_NN, double* NEO_LN,
			      double* NEO_YN, int* PHASE, int* PHASE_MID, int* PHASE_COPT, double* PHASE_CORE, int* PHASE_INTN, int* PHASE_INTU, int* PHASE_STAGE5, int* PHASE_OLD,
			      int* PHASE_FREQ, int* PHASE_LIN, int* PHASE_HIGH, int* PHASE_RATS, double* PHASE_SCALE, double* PHASE_CHIR, int* RESTART,
			      double* TSTART, double* TEND, double* DT); 

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
  int    NEO_INTP;       // If != 0 then interpolate pFiles in NEOCLASSICAL
  int	 NEO_INTF;       // If != 0 then interpolate fFiles in NEOCLASSICAL
  int	 NEO_IMPURITY;   // If != 0 then include single ion impurity species in NEOCLASSICAL
  int	 NEO_NEUTRAL;    // If != 0  then include majority ion neutrals in NEOCLASSICAL
  int	 NEO_FREQ;       // Set natural frequency type (see NEOCLASSICAL/Inputs/Neoclassical.nml) in NEOCLASSICAL
  int    NEO_NTYPE;      // If 0/1 then neutral density distribution exponential/Lorentzian in NEOCLASSICAL
  double NEO_NN;         // Set neutral density at boundary in NEOCLASSICAL
  double NEO_LN;         // Set neutral density scalelength in NEOCLASSICAL
  double NEO_YN;         // Set neutral poloidal asymmetry parameter in NEOCLASSICAL
  int    PHASE;          // If != 0 then call PHASE
  int	 PHASE_MID;      // If != 0 then middle RMP coils included in PHASE calculation
  int    PHASE_COPT;     // If != 0 then optimize coil currents
  double PHASE_CORE;     // Core drive minimization factor
  int	 PHASE_INTF;     // If != 0 then interpolate fFiles in PHASE
  int	 PHASE_INTN;     // If != 0 then interpolate nFiles in PHASE
  int	 PHASE_INTU;     // If != 0 then interpolate uFiles/mFiles/lFiles in PHASE
  int	 PHASE_STAGE5;   // If != 0 then Stage5 PHASE calculation enabled
  int	 PHASE_OLD;      // If != 0 then restart PHASE calculations from previous run
  int    PHASE_FREQ;     // If != 0 then use island width dependent natural frequency in PHASE (overrides NEO_FREQ)
  int    PHASE_LIN;      // If != 0 then perform purely linear PHASE calculation 
  int    PHASE_HIGH;     // If != 0 then perfprm higher order transport calculation in PHASE
  int    PHASE_RATS;     // If != 0 then perform linear-only interpolation of uFiles/mFiles/lFiles in PHASE
  double PHASE_SCALE;    // GPEC scalefactor in PHASE
  double PHASE_CHIR;     // Maximum Chirikov parameter for vacuum islands in PHASE
  int	 RESTART;        // If != 0 then delete all previous IslandDynamics data
  double TSTART;         // Simulation experimental start time (ms)
  double TEND;           // Simulation experimental end time (ms)
  double DT;             // Simulation experimental time step (ms)

  NameListRead (&FLUX_NTOR, &FLUX_MMIN, &FLUX_MMAX,
		&NEO_INTP, &NEO_INTF, &NEO_IMPURITY, &NEO_NEUTRAL, &NEO_FREQ, &NEO_NTYPE, &NEO_NN, &NEO_LN, &NEO_YN,
		&PHASE, &PHASE_MID, &PHASE_COPT, &PHASE_CORE, &PHASE_INTN, &PHASE_INTU, &PHASE_STAGE5, &PHASE_OLD, &PHASE_FREQ, &PHASE_LIN,
		&PHASE_HIGH, &PHASE_RATS, &PHASE_SCALE, &PHASE_CHIR, &RESTART, &TSTART, &TEND, &DT);

  PHASE_INTF = NEO_INTF;

  printf ("Reading Inputs/Island.nml:\n");
  printf ("FLUX_NTOR = %2d  FLUX_MMIN = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  printf ("NEO_INTP = %2d  NEO_INTF = %2d  NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  printf ("PHASE = %2d  PHASE_MID = %2d  PHASE_COPT = %2d  PHASE_CORE = %11.4e  PHASE_INTN = %2d  PHASE_INTU = %2d  PHASE_STAGE5 = %2d  PHASE_OLD = %2d  PHASE_FREQ = %2d  PHASE_LIN = %2d  PHASE_HIGH = %2d  PHASE_RATS = %2d  PHASE_SCALE = %11.4e  PHASE_CHIR = %11.4e\n",
	  PHASE, PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, PHASE_HIGH, PHASE_RATS, PHASE_SCALE, PHASE_CHIR);
  printf ("RESTART = %2d  TSTART = %11.4e  TEND = %11.4e  DT = %11.4e  IRMP = %11.4e\n",
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
  printf ("FLUX_NTOR = %2d  FLUX_MMIN = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  printf ("NEO_INTP = %2d  NEO_INTF = %2d  NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	  NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  printf ("PHASE = %2d  PHASE_MID = %2d  PHASE_COPT = %2d  PHASE_CORE = %11.4e  PHASE_INTN = %2d  PHASE_INTU = %2d  PHASE_STAGE5 = %2d  PHASE_OLD = %2d  PHASE_FREQ = %2d  PHASE_LIN = %2d  PHASE_HIGH = %2d  PHASE_RATS = %2d  PHASE_SCALE = %11.4e  PHASE_CHIR = %11.4e\n",
	  PHASE, PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, PHASE_HIGH, PHASE_RATS, PHASE_SCALE, PHASE_CHIR);
  printf ("RESTART = %2d  TSTART = %11.4e  TEND = %11.4e  DT = %11.4e  IRMP = %11.4e\n",
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
  fprintf (monitor, "FLUX_NTOR = %2d  FLUX_MMIN  = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  fprintf (monitor, "NEO_INTP = %2d  NEO_INTF = %2d  NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	   NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  fprintf (monitor, "PHASE = %2d  PHASE_MID = %2d  PHASE_COPT = %2d  PHASE_CORE = %11.4e  PHASE_INTN = %2d  PHASE_INTU = %2d  PHASE_STAGE5 = %2d  PHASE_OLD = %2d  PHASE_FREQ = %2d  PHASE_LIN = %2d  PHASE_HIGH = %2d  PHASE_RATS = %2d  PHASE_SCALE = %11.4e  PHASE_CHIR = %11.4e\n",
	   PHASE, PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, PHASE_HIGH, PHASE_RATS, PHASE_SCALE, PHASE_CHIR);
  fprintf (monitor, "RESTART = %2d  TSTART = %11.4e  TEND = %11.4e  DT = %11.4e  IRMP = %11.4e\n",
	   RESTART, TSTART, TEND, DT, IRMP);
  fclose (monitor);

  FILE* namelist = fopen ("Inputs/InputParameters", "w");
  fprintf (namelist, "Input parameters (from Inputs/Island.nml):\n");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "FLUX_NTOR = %2d  FLUX_MMIN = %2d  FLUX_MMAX = %2d\n",
	  FLUX_NTOR, FLUX_MMIN, FLUX_MMAX);
  fprintf (namelist, "NEO_INTP = %2d  NEO_INTF = %2d  NEO_IMPURITY = %2d  NEO_NEUTRAL = %2d  NEO_FREQ = %2d  NEO_NTYPE = %2d  NEO_NN = %11.4e  NEO_LN = %11.4e  NEO_YN = %11.4e\n",
	   NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN);
  fprintf (namelist, "PHASE = %2d  PHASE_MID = %2d  PHASE_COPT = %2d  PHASE_CORE = %11.4e  PHASE_INTN = %2d  PHASE_INTU = %2d  PHASE_STAGE5 = %2d  PHASE_OLD = %2d  PHASE_FREQ = %2d  PHASE_LIN = %2d  PHASE_HIGH = %2d  PHASE_RATS = %2d  PHASE_SCALE = %11.4e  PHASE_CHIR = %11.4e\n",
	   PHASE, PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, PHASE_HIGH, PHASE_RATS, PHASE_SCALE, PHASE_CHIR);
  fprintf (namelist, "RESTART = %2d  TSTART = %11.4e  TEND = %11.4e  DT = %11.4e  IRMP = %11.4e\n",
	   RESTART, TSTART, TEND, DT, IRMP);
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
	sprintf (PHASEstring, "cd ../Phase; ./phase -m %d -C %d -f %d -n %d -u %d -s %2d -o %2d -F %2d -l %2d -t %16.9e -T %16.9e -S %16.9e -c %11.4e -i %11.4e -H %2d -r %2d -D %16.9e",
		 PHASE_MID, PHASE_COPT, PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, Time, Time + DT, PHASE_SCALE, PHASE_CHIR, IRMP, PHASE_HIGH, PHASE_RATS, PHASE_CORE);
      else
	sprintf (PHASEstring, "cd ../Phase; ./phase -m %d -C %d -f %d -n %d -u %d -s %2d -o %2d -F %2d -l %2d -t %16.9e -T %16.9e -S %16.9e -c %11.4e -i %11.4e -H %2d -r %2d -D %16.9e",
		 PHASE_MID, PHASE_COPT, PHASE_INTF, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, lock, PHASE_FREQ, PHASE_LIN, Time, Time + DT, PHASE_SCALE, PHASE_CHIR, IRMP, PHASE_HIGH, PHASE_RATS, PHASE_CORE);
      
      lock = 1;
      
      // Call program FLUX
      if (NEO_INTF == 0)
	{
	  printf ("Executing:: %s\n", FLUXstring);
	  if (CallSystem (FLUXstring) != 0)
	    exit (1);
	}

      // Call program NEOCLASSICAL
      if (PHASE_INTN == 0)
	{
	  printf ("Executing:: %s\n", NEOstring);
	  if (CallSystem (NEOstring) != 0)
	    exit (1);
	}

      // Call program PHASE
      if (PHASE)
	{
	  printf ("Executing:: %s\n", PHASEstring);
	  if (CallSystem (PHASEstring) != 0)
	    exit (1);
	}

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
