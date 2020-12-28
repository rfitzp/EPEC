// main.cpp

// ###############################
// Main function for program Phase
// See Phase.h
// ###############################

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Phase.h"

int main (int argc, char** argv)
{
  // .....................
  // Print welcome message
  // .....................
  printf ("\n#############\nProgram PHASE\n#############\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  FILE* monitor = fopen ("../IslandDynamics/Outputs/monitor.txt", "a");
  fprintf (monitor, "\n#############\nProgram PHASE\n#############\n");
  fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  fclose (monitor);
  
  // ........................
  // Get command line options
  // ........................
  int c;
  char* tvalue = NULL; char* svalue = NULL; char* fvalue = NULL; 
  char* nvalue = NULL; char* uvalue = NULL; char* ovalue = NULL;
  char* Fvalue = NULL; char* Svalue = NULL; char* lvalue = NULL;
  char* mvalue = NULL; char* Tvalue = NULL;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "f:hn:o:s:t:u:F:S:l:m:T:")) != -1)
    switch (c)
      {
      case 'f':
	fvalue = optarg;
 	break;
      case 'h':
 	printf ("Options:\n");
	printf ("-f INTF   - set interpolation flag INTF\n");
	printf ("-h        - list options\n");
	printf ("-l LIN    - set linear flag LIN\n");
	printf ("-m MID    - set mFile flag MID\n");
	printf ("-n INTN   - set interpolation flag INTN\n");
	printf ("-o OLD    - set old calculation flag OLD\n");
	printf ("-s STAGE5 - set Stage5 flag STAGE5\n");
	printf ("-t TSTART - set simulation start time (s) to TSTART\n");
	printf ("-T TEND   - set simulation end time (s) to TEND\n");
	printf ("-u INTU   - set interpolation flag INTU\n");
	printf ("-F FREQ   - set frequency flag FREQ\n");
	printf ("-S SCALE  - set GPEC scalefactor SCALE\n");
	exit (0);
      case 'n':
	nvalue = optarg;
 	break;
      case 'o':
	ovalue = optarg;
 	break;
      case 's':
	svalue = optarg;
 	break;
      case 't':
	tvalue = optarg;
 	break;
      case 'T':
	Tvalue = optarg;
 	break;
      case 'u':
	uvalue = optarg;
 	break;
      case 'l':
	lvalue = optarg;
 	break;
      case 'm':
	mvalue = optarg;
 	break;
      case 'F':
	Fvalue = optarg;
 	break;
      case 'S':
	Svalue = optarg;
 	break;
      case '?':
	if (optopt == 't' || optopt == 's' || optopt == 'f' || optopt == 'o'|| optopt == 'n' || optopt == 'u' || optopt == 'F' || optopt == 'S' || optopt == 'l' || optopt == 'm' || optopt == 'T')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _STAGE5 = -1;    int   _INTF = -1; int    _OLD    = -1;    int    _FREQ = -1;
  int    _INTN   = -1;    int   _INTU = -1; double _TSTART = -1.e6; float  __TSTART;
  double _SCALE  = -1.e6; float __SCALE;    int    _LIN    = -1.;   int    _MID = -1;
  double _TEND   = -1.e6; float __TEND;
  
  if (svalue != NULL)
    _STAGE5 = atoi (svalue);
  if (fvalue != NULL)
    _INTF= atoi (fvalue);
  if (nvalue != NULL)
    _INTN= atoi (nvalue);
  if (uvalue != NULL)
    _INTU= atoi (uvalue);
  if (ovalue != NULL)
     _OLD = atoi (ovalue);
  if (lvalue != NULL)
    _LIN = atoi (lvalue);
  if (mvalue != NULL)
    _MID = atoi (mvalue);
  if (Fvalue != NULL)
     _FREQ = atoi (Fvalue);
  if (tvalue != NULL)
    {
      __TSTART = atof (tvalue);
      _TSTART = double (__TSTART);
    }
   if (Tvalue != NULL)
    {
      __TEND = atof (Tvalue);
      _TEND  = double (__TEND);
    }
   if (Svalue != NULL)
    {
      __SCALE = atof (Svalue);
      _SCALE  = double (__SCALE);
    }
    
  // ..................
  // Call program PHASE
  // ..................
  Phase phase;
  phase.Solve (_STAGE5, _INTF, _INTN, _INTU, _OLD, _FREQ, _LIN, _MID, _TSTART, _TEND, _SCALE);

  return 0;
}
