// main.cpp

// ######################################
// Main function for program Neoclassical
// See Neoclassical.h
// ######################################

#include "Neoclassical.h"

int main (int argc, char** argv)
{
  // Print welcome message
  printf ("\n####################\nProgram NEOCLASSICAL\n####################\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

   // Get options
  int c;
  char* nvalue = NULL; char* Ivalue = NULL; char* fvalue = NULL; char* tvalue = NULL;
  char* yvalue = NULL; char* ivalue = NULL; 
  opterr = 0;
  
  while ((c = getopt (argc, argv, "f:hin:t:y:I:")) != -1)
    switch (c)
      {
      case 'f':
	fvalue = optarg;
	break;
      case 'h':
	printf ("Options:\n");
	printf ("-f FREQ     - set frequency flag\n");
	printf ("-h          - list options\n");
	printf ("-i INTP     - set interpolation flag\n");
	printf ("-n NEUTRAL  - set neutral flag\n");
	printf ("-t TIME     - set experimental time\n");
	printf ("-y YN       - set neutral peaking factor\n");
	printf ("-I IMPURITY - set impurity flag\n");
	exit (0);
      case 'i':
	ivalue = optarg;
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
	if (optopt == 'n' || optopt == 'I' || optopt == 'f' || optopt == 't' || optopt == 'y' || optopt == 'i')
	  printf ("Opetion = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _NEUTRAL = -1; int _IMPURITY = -1; int _FREQ = -999999; int _INTP = -1;
  double _TIME = 0.; double _YN = -1.;

  if (nvalue != NULL)
    _NEUTRAL = atoi (nvalue);
  if (Ivalue != NULL)
    _IMPURITY = atoi (Ivalue);
  if (fvalue != NULL)
    _FREQ = atoi (fvalue);
  if (ivalue != NULL)
    _INTP = atoi (ivalue);
  if (yvalue != NULL)
    _YN = double (atof (yvalue));
  if (tvalue != NULL)
    _TIME = double (atof (tvalue));
 
  // Call program
  Neoclassical neoclassical;
  neoclassical.Solve (_NEUTRAL, _IMPURITY, _FREQ, _INTP, _YN, _TIME);

  return 0;
}
