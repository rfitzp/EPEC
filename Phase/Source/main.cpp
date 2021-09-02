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

  // ..................
  // Call program PHASE
  // ..................
  Phase phase;
  clock_t begin = clock ();
  phase.Solve ();
  clock_t end = clock ();
  double time_spent = double (end - begin) /double (CLOCKS_PER_SEC);

  // ..................
  // Print exit message
  // ..................
  printf ("*************************************************************\n");
  printf ("PROGRAM PHASE:: Normal termination: Wall time = %11.4e s\n", time_spent);
  printf ("*************************************************************\n");

  return 0;
}
