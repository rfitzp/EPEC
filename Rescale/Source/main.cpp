// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Rescale
// See Rescale.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Rescale.h"

int main (int argc, char** argv)
{
  // ...............
  // Welcome message
  // ...............
  printf ("\n###############\nProgram RESCALE\n###############\n");
  printf ("Version: %1d.%1d\n", VERSION_MAJOR, VERSION_MINOR);
  
  // ....................
  // Call program RESCALE
  // ....................
  Rescale rescale;
  clock_t begin = clock ();
  rescale.RescaleEquilibrium ();

  clock_t end = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  // ..................
  // Print exit message
  // ..................
  printf ("***************************************************************\n");
  printf ("PROGRAM RESCALE:: Normal termination: Wall time = %11.4e s\n", time_spent);
  printf ("***************************************************************\n");

  return 0;
}
