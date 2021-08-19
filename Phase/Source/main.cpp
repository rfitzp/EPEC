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

  // ........................
  // Get command line options
  // ........................
  int c;
  char* tvalue = NULL; char* svalue = NULL; char* fvalue = NULL; 
  char* nvalue = NULL; char* uvalue = NULL; char* ovalue = NULL;
  char* Fvalue = NULL; char* Svalue = NULL; char* lvalue = NULL;
  char* mvalue = NULL; char* Tvalue = NULL; char* cvalue = NULL;
  char* ivalue = NULL; char* rvalue = NULL; char* Hvalue = NULL;
  char* Cvalue = NULL; char* Dvalue = NULL; char* Nvalue = NULL;
  char* Gvalue = NULL; char* Xvalue = NULL; char* Bvalue = NULL;
  char* Vvalue = NULL; char* Pvalue = NULL;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "f:hn:o:s:t:u:S:l:m:T:c:i:r:H:C:D:F:N:G:X:B:V:P:")) != -1)
    switch (c)
      {
      case 'f':
	fvalue = optarg;
 	break;
      case 'h':
 	printf ("Options:\n");
	printf ("-c CHIR   - set maximum vacuum Chirkov parameter to CHIR\n");
	printf ("-f INTF   - set fFile interpolation enabling flag to INTF\n");
	printf ("-h        - list options\n");
	printf ("-i IRMP   - set RMP current to IRMP");
	printf ("-l LIN    - set linear calculation enabling flag to LIN\n");
	printf ("-m MID    - set mFile enabling flag to MID\n");
	printf ("-n INTN   - set nFile interpolation enabling flag to INTN\n");
	printf ("-o OLD    - set old calculation flag to OLD\n");
	printf ("-r RATS   - set uFile/mFile/lFile linear-only interpolation flag to RATS\n");
	printf ("-s STAGE5 - set Stage5 caculation enabling flag to STAGE5\n");
	printf ("-t TSTART - set simulation start time (s) to TSTART\n");
	printf ("-u INTU   - set uFile/mFile/lFile interpolation enabling flag to INTU\n");
	printf ("-B BOOT   - set bootstrap flag to BOOT\n");
	printf ("-C COPT   - set coil current optimization flag to COPT\n");
	printf ("-D CORE   - set core minimization factor to CORE\n");
	printf ("-F FREQ   - set natural frequency selection flag to FREQ\n");
	printf ("-G FFAC   - set natural frequency selection parameter to FFAC\n");
	printf ("-H HIGH   - set higher order transport calculation enabling flag to HIGH\n");
	printf ("-P POLZ   - set polarization flag to POLZ\n");
	printf ("-N NATS   - set linear only nFiles interpolation flag to NATS\n");
	printf ("-S SCALE  - set GPEC scalefactor to SCALE\n");	
	printf ("-T TEND   - set simulation end time (s) to TEND\n");
	printf ("-V CURV   - set curvature flag to CURV\n");
	printf ("-X CXD    - set charge exchange damping flag to CXD\n");
	exit (0);
      case 'n':
	nvalue = optarg;
 	break;
      case 'o':
	ovalue = optarg;
 	break;
      case 'i':
	ivalue = optarg;
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
      case 'F':
	Fvalue = optarg;
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
      case 'S':
	Svalue = optarg;
 	break;
      case 'c':
	cvalue = optarg;
 	break;
      case 'r':
	rvalue = optarg;
 	break;
      case 'H':
	Hvalue = optarg;
 	break;
      case 'C':
	Cvalue = optarg;
 	break;
      case 'D':
	Dvalue = optarg;
 	break;
      case 'N':
	Nvalue = optarg;
 	break;
      case 'G':
	Gvalue = optarg;
 	break;
      case 'X':
	Xvalue = optarg;
 	break;
      case 'B':
	Bvalue = optarg;
 	break;
      case 'V':
	Vvalue = optarg;
 	break;
      case 'P':
	Pvalue = optarg;
 	break;
      case '?':
	if (optopt == 't' || optopt == 's' || optopt == 'f' || optopt == 'o'|| optopt == 'n' || optopt == 'u'
	    || optopt == 'S' || optopt == 'l' || optopt == 'm' || optopt == 'T' || optopt == 'c'
	    || optopt == 'i' || optopt == 'r' || optopt == 'H' || optopt == 'C' || optopt == 'D' || optopt == 'F'|| optopt == 'N'
	    || optopt == 'X' || optopt == 'B' || optopt == 'V' || optopt == 'P')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _STAGE5 = -1;    int   _INTF = -1;     int    _OLD    = -1;    int   _FREQ = -1;    
  int    _INTN   = -1;    int   _INTU = -1;     double _TSTART = -1.e6; float __TSTART;
  double _SCALE  = -1.e6; float __SCALE;        int    _LIN    = -1.;   int   _MID  = -1;
  double _TEND   = -1.e6; float __TEND;         double _CHIR   = -1.e6; float __CHIR;
  double _IRMP   = -1.e6; float __IRMP;         int    _HIGH   = -1.;   int   _RATS = -1.;
  int    _COPT   = -1;    double _CORE = -1.e6; float __CORE;           int   _NATS = -1.;
  double _FFAC   = -1.e6; float __FFAC;         int    _CXD    = -1;    int   _BOOT = -1;
  int    _CURV   = 1;     int    _POLZ = -1.;
  
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
  if (Hvalue != NULL)
     _HIGH = atoi (Hvalue);
  if (rvalue != NULL)
     _RATS = atoi (rvalue);
  if (Nvalue != NULL)
     _NATS = atoi (Nvalue);
  if (Cvalue != NULL)
     _COPT = atoi (Cvalue);
  if (Fvalue != NULL)
     _FREQ = atoi (Fvalue);
  if (Xvalue != NULL)
    _CXD = atoi (Xvalue);
  if (Bvalue != NULL)
     _BOOT = atoi (Bvalue);
  if (Vvalue != NULL)
     _CURV = atoi (Vvalue);
  if (Pvalue != NULL)
     _POLZ = atoi (Pvalue);
  if (tvalue != NULL)
    {
      __TSTART = atof (tvalue);
      _TSTART  = double (__TSTART);
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
   if (cvalue != NULL)
     {
      __CHIR = atof (cvalue);
      _CHIR  = double (__CHIR);
     }
   if (ivalue != NULL)
     {
      __IRMP = atof (ivalue);
      _IRMP  = double (__IRMP);
     }
   if (Dvalue != NULL)
     {
      __CORE = atof (Dvalue);
      _CORE  = double (__CORE);
     }
   if (Gvalue != NULL)
     {
      __FFAC = atof (Gvalue);
      _FFAC  = double (__FFAC);
     }
    
  // ..................
  // Call program PHASE
  // ..................
  Phase phase;
  phase.Solve (_STAGE5, _INTF, _INTN, _INTU, _NATS, _OLD, _LIN, _MID, _COPT, _TSTART, _TEND,
	       _SCALE, _CHIR, _IRMP, _HIGH, _RATS, _CORE, _FREQ, _FFAC, _CXD, _BOOT, _CURV, _POLZ);

  return 0;
}
