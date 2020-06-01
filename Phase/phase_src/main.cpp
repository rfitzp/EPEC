// main.cpp

// ###############################
// Main function for program Phase
// ###############################

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Phase.h"

int main (int argc, char** argv)
{
  // Print welcome message
  printf ("\n#############\nProgram PHASE\n#############\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // Get options
  int c;
  char* tvalue = NULL; char* svalue = NULL; char* ivalue = NULL; char* ovalue = NULL;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "hi:o:s:t:")) != -1)
    switch (c)
      {
      case 'h':
 	printf ("Options:\n");
	printf ("-h        - list options\n");
	printf ("-i INTP   - set interpolation flag\n");
	printf ("-o OLD    - set old calculation flag\n");
	printf ("-s STAGE2 - set Stage2 flag\n");
	printf ("-t TIME   - set time\n");
	exit (0);
      case 'i':
	ivalue = optarg;
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
      case '?':
	if (optopt == 't' || optopt == 's' || optopt == 'i' || optopt == 'o')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _STAGE2 = -1; int _INTP = -1; int _OLD = -1;
  float  __TIME;  double _TIME = -1.e6;
  
  if (svalue != NULL)
    _STAGE2 = atoi (svalue);
  if (ivalue != NULL)
    _INTP = atoi (ivalue);
  if (ovalue != NULL)
     _OLD = atoi (ovalue);
  if (tvalue != NULL)
    {
      __TIME = atof (tvalue);
      _TIME  = double (__TIME);
    }

  // Call program
  Phase phase;
  phase.Solve (_STAGE2, _INTP, _OLD, _TIME);

  return 0;
}
