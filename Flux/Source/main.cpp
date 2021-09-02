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

int pRhs6 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs6 (r, y, dydr, NULL);

  return status;
}

int pRhs7 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs7 (r, y, dydr, NULL);

  return status;
}

int pRhs8 (double r, const double y[], double dydr[], void* params)
{
  Flux flux = *(Flux*) params;

  int status = flux.Rhs8 (r, y, dydr, NULL);

  return status;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main (int argc, char** argv)
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n############\nProgram FLUX\n############\n");
  printf ("Version: %1d.%1d\n", VERSION_MAJOR, VERSION_MINOR);
  
  // .................
  // Call program FLUX
  // .................
  Flux flux;
  clock_t begin = clock ();
  flux.Solve ();

  clock_t end = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  // ..................
  // Print exit message
  // ..................
  printf ("*************************************************************\n");
  printf ("PROGRAM FLUX:: Normal termination: Wall time = %11.4e s\n", time_spent);
  printf ("*************************************************************\n");

  return 0;
}
