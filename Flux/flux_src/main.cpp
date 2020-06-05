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
  // ...............
  // Welcome message
  // ...............
  printf ("\n############\nProgram FLUX\n############\n");
  printf ("Version: %1d.%1d\n\n", VERSION_MAJOR, VERSION_MINOR);

  // ........................
  // Get command line options
  // ........................
  int c;
  char* nvalue = NULL; char* mvalue = NULL; char* Mvalue = NULL;
  char* tvalue = NULL; char* gvalue = NULL;
  opterr = 0;
  
  while ((c = getopt (argc, argv, "hg:n:m:t:M:")) != -1)
    switch (c)
      {
      case 'h':
	printf ("Options:\n");
	printf ("-h      - list options\n");
	printf ("-g INTG - set interpolation flag INTG\n");
	printf ("-n NTOR - set toroidal mode number to NTOR\n");
	printf ("-m MMIN - set minumum mode number to MMIN\n");
	printf ("-t TIME - set experimental time to TIME\n");
	printf ("-M MMAX - set maximum mode number to MMAX\n");
	exit (0);
     case 'g':
	gvalue = optarg;
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
	if (optopt == 'n' || optopt == 'm' || optopt == 'M' || optopt == 't' || optopt == 'g')
	  printf ("Option = %c requires an argument\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option '-%c'\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'\n", optopt);
	return 1;
      default:
	abort ();
      }

  int    _NTOR = -1, _MMIN = -1, _MMAX = -1, _INTG = -1;
  double _TIME = 0.;
  
  if (nvalue != NULL)
    _NTOR = atoi (nvalue);
  if (mvalue != NULL)
    _MMIN = atoi (mvalue);
  if (Mvalue != NULL)
    _MMAX = atoi (Mvalue);
  if (tvalue != NULL)
    _TIME = double (atof (tvalue));
  if (gvalue != NULL)
    _INTG = atoi (gvalue);

  // .................
  // Call program FLUX
  // .................
  Flux flux;
  flux.Solve (_INTG, _NTOR, _MMIN, _MMAX, _TIME);
 
  return 0;
}
