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

  // .........................
  // Call program NEOCLASSICAL
  // .........................
  Neoclassical neoclassical;
  clock_t begin = clock ();
  neoclassical.Solve ();
  clock_t end = clock ();
  double time_spent = double (end - begin) /double (CLOCKS_PER_SEC);

  // ..................
  // Print exit message
  // ..................
  printf ("********************************************************************\n");
  printf ("PROGRAM NEOCLASSICAL:: Normal termination: Wall time = %11.4e s\n", time_spent);
  printf ("********************************************************************\n");

  return 0;
}
