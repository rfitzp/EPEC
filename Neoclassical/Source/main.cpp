// main.cpp

// ######################################
// Main function for program Neoclassical
// See Neoclassical.h
// ######################################

#include "Neoclassical.h"

int main (int argc, char** argv)
{
  // .....................
  // Print welcome message
  // .....................
  printf ("\n####################\nProgram NEOCLASSICAL\n####################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  FILE* monitor = fopen ("../IslandDynamics/Outputs/monitor.txt", "a");
  fprintf (monitor, "\n####################\nProgram NEOCLASSICAL\n####################\n");
  fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
  fclose (monitor);
			   
  // ........................
  // Get command line options
  // ........................
  int c;
  char* nvalue = NULL; char* Ivalue = NULL; char* fvalue = NULL; char* tvalue = NULL;
  char* yvalue = NULL; char* pvalue = NULL; char* evalue = NULL; 
  opterr = 0;
  
  while ((c = getopt (argc, argv, "e:f:hp:n:t:y:I:")) != -1)
    switch (c)
      {
      case 'e':
	evalue = optarg;
	break;
      case 'f':
	fvalue = optarg;
	break;
      case 'h':
	printf ("Options:\n");
	printf ("-e INTF     - set interpolation flag INTF\n");
	printf ("-f FREQ     - set frequency flag FREQ\n");
	printf ("-h          - list options\n");
	printf ("-n NEUTRAL  - set neutral flag NEUTRAL\n");
	printf ("-p INTP     - set interpolation flag INTP\n");
	printf ("-t TIME     - set experimental time to TIME\n");
	printf ("-y YN       - set neutral peaking factor YN\n");
	printf ("-I IMPURITY - set impurity flag IMPURITY\n");
	exit (0);
      case 'p':
	pvalue = optarg;
	break;
      case 'n':
	nvalue = optarg;
 	break;
      case 't':
	tvalue = optarg;
	break;
      case 'y':
	yvalue = optarg;
	break;
      case 'I':
	Ivalue = optarg;
	break;
      case '?':
	if (optopt == 'n' || optopt == 'I' || optopt == 'f' || optopt == 't' || optopt == 'y' || optopt == 'p' || optopt == 'e')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _NEUTRAL = -1;  int    _IMPURITY = -1; int _FREQ = -999999; int _INTP = -1;
  double _TIME    = -1.; double _YN       = -1.;int _INTF = -1;

  if (nvalue != NULL)
    _NEUTRAL = atoi (nvalue);
  if (Ivalue != NULL)
    _IMPURITY = atoi (Ivalue);
  if (fvalue != NULL)
    _FREQ = atoi (fvalue);
  if (pvalue != NULL)
    _INTP = atoi (pvalue);
  if (evalue != NULL)
    _INTF = atoi (evalue);
  if (yvalue != NULL)
    _YN = double (atof (yvalue));
  if (tvalue != NULL)
    _TIME = double (atof (tvalue));

  // .........................
  // Call program NEOCLASSICAL
  // .........................
  Neoclassical neoclassical;
  neoclassical.Solve (_NEUTRAL, _IMPURITY, _FREQ, _INTP, _INTF, _YN, _TIME);

  return 0;
}
