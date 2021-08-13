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
	   
  // ........................
  // Get command line options
  // ........................
  int c;
  char* nvalue = NULL; char* Ivalue = NULL; char* fvalue = NULL; char* tvalue = NULL;
  char* yvalue = NULL; char* pvalue = NULL; char* evalue = NULL; char* lvalue = NULL;
  char* Nvalue = NULL; char* Tvalue = NULL; char* cvalue = NULL; char* Cvalue = NULL; 
  int   _OMFIT = 0;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "c:e:f:hp:n:t:y:I:l:N:T:C:o")) != -1)
    switch (c)
      {
      case 'c':
	cvalue = optarg;
	break;
      case 'e':
	evalue = optarg;
	break;
      case 'f':
	fvalue = optarg;
	break;
      case 'h':
	printf ("Options:\n");
	printf ("-c INTC     - set interpolation flag to INTC\n");
	printf ("-e INTF     - set interpolation flag to INTF\n");
	printf ("-f EXB      - set ExB frequency flag to EXB\n");
	printf ("-h          - list options\n");
	printf ("-l LN       - set neutral density distribution scalelength to LN\n");
	printf ("-n NEUTRAL  - set neutral flag to NEUTRAL\n");
	printf ("-o          - flag to select OMFIT mode\n");
	printf ("-p INTP     - set interpolation flag to INTP\n");
	printf ("-t TIME     - set experimental time to TIME\n");
	printf ("-y YN       - set neutral poloidal peaking factor to YN\n");
	printf ("-C CATS     - set cFile interpolation flag to CATS\n");
	printf ("-I IMPURITY - set impurity flag to IMPURITY\n");
	printf ("-N NN       - set boundary neutral density to NN\n");
	printf ("-T NTYPE    - set neutral density distribution type to NTYPE\n");
	exit (0);
      case 'o':
	_OMFIT = 1;
      case 'p':
	pvalue = optarg;
	break;
      case 'l':
	lvalue = optarg;
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
      case 'N':
	Nvalue = optarg;
	break;
      case 'T':
	Tvalue = optarg;
	break;
      case 'C':
	Cvalue = optarg;
	break;
      case '?':
	if (optopt == 'n' || optopt == 'I' || optopt == 'f' || optopt == 't' || optopt == 'y' || optopt == 'p' || optopt == 'e'
	    || optopt == 'l' || optopt == 'N' || optopt == 'T' || optopt == 'c' || optopt == 'C')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _NEUTRAL = -1;  int    _IMPURITY = -1;  int _EXB   = -999999; int _INTP = -1;
  double _TIME    = -1.; double _YN       = -1.; int _INTF  = -1;      int _INTC = -1;
  double _LN      = -1.; double _NN       = -1.; int _NTYPE = -1;      int _CATS = -1;

  if (nvalue != NULL)
    _NEUTRAL = atoi (nvalue);
  if (Ivalue != NULL)
    _IMPURITY = atoi (Ivalue);
  if (fvalue != NULL)
    _EXB = atoi (fvalue);
  if (pvalue != NULL)
    _INTP = atoi (pvalue);
  if (evalue != NULL)
    _INTF = atoi (evalue);
  if (yvalue != NULL)
    _YN = double (atof (yvalue));
  if (tvalue != NULL)
    _TIME = double (atof (tvalue));
  if (lvalue != NULL)
    _LN = double (atof (lvalue));
  if (Nvalue != NULL)
    _NN = double (atof (Nvalue));
  if (Tvalue != NULL)
    _NTYPE = atoi (Tvalue);
  if (cvalue != NULL)
    _INTC = atoi (cvalue);
  if (Cvalue != NULL)
    _CATS = atoi (Cvalue);

  if (_OMFIT)
    printf ("OMFIT mode\n\n");
  else
    {
      printf ("Normal mode\n\n");
     
      FILE* monitor = fopen ("../IslandDynamics/Outputs/monitor.txt", "a");
      fprintf (monitor, "\n####################\nProgram NEOCLASSICAL\n####################\n");
      fprintf (monitor, "Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);
      fprintf (monitor, "Normal mode\n\n");
      fclose (monitor);
    }
  fflush (stdout);

  // Call program NEOCLASSICAL
  // .........................
  Neoclassical neoclassical;
  neoclassical.Solve (_NEUTRAL, _IMPURITY, _EXB, _INTP, _INTF, _INTC, _NTYPE, _NN, _LN, _YN, _TIME, _CATS, _OMFIT);

  return 0;
}
