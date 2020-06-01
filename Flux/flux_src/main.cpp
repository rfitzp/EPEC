// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Flux
// See Flux.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Flux.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Pointers to right-hand side functions for passing to gsl adaptive integration routines
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int pRhs1 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs1 (r, y, dydr, NULL);

  return status;
}

int pRhs2 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs2 (r, y, dydr, NULL);

  return status;
}

int pRhs3 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs3 (r, y, dydr, NULL);

  return status;
}

int pRhs4 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs4 (r, y, dydr, NULL);

  return status;
}

int pRhs5 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs5 (r, y, dydr, NULL);

  return status;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main (int argc, char** argv)
{
  // Welcome message
  printf ("\n############\nProgram FLUX\n############\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // Get options
  int c;
  char* nvalue = NULL; char* mvalue = NULL; char* Mvalue = NULL;
  char* tvalue = NULL; char* ivalue = NULL;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "hi:n:m:t:M:")) != -1)
    switch (c)
      {
      case 'h':
	printf ("Options:\n");
	printf ("-h      - list options\n");
	printf ("-i INTP - set interpolation flag\n");
	printf ("-n NTOR - set toroidal mode number\n");
	printf ("-m MMIN - set minumum mode number\n");
	printf ("-t TIME - set experimental time\n");
	printf ("-M MMAX - set maximum mode number\n");
	exit (0);
     case 'i':
	ivalue = optarg;
	break;
     case 'n':
	nvalue = optarg;
 	break;
      case 'm':
	mvalue = optarg;
	break;
      case 't':
	tvalue = optarg;
	break;
      case 'M':
	Mvalue = optarg;
	break;
      case '?':
	if (optopt == 'n' || optopt == 'm' || optopt == 'M' || optopt == 't' || optopt == 'i')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _NTOR = -1, _MMIN = -1, _MMAX = -1, _INTP = -1;
  double _TIME = 0.;
  
  if (nvalue != NULL)
    _NTOR = atoi (nvalue);
  if (mvalue != NULL)
    _MMIN = atoi (mvalue);
  if (Mvalue != NULL)
    _MMAX = atoi (Mvalue);
  if (tvalue != NULL)
    _TIME = double (atof (tvalue));
  if (ivalue != NULL)
    _INTP = atoi (ivalue);

  // Call program
  Flux flux;
  flux.Solve (_INTP, _NTOR, _MMIN, _MMAX, _TIME);
 
  return 0;
}
