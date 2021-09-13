// Rescale.h

// ##############################################
// Class to rescale equilibrium gFile and pFile

// Rescaling such that q_95 = q_95_new, where
// q95_new is read from Inputs/Rescale.nml, but
// vacuum toroidal magnetic field remains constant

// Initial gFile in Inputs/gFile
// Rescaled gFile in Outputs/gFile

// Initial pFile in Inputs/pFile
// Rescaled pFile in Outputs/pFile

// ##############################################

#ifndef RESCALE
#define RESCALE

#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define MAXFILENAMELENGTH 500

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include "Field.h"

// gFile rescaling function
extern "C" void gFileRescale (double* q95_old, double* q95_new, double* a1);

// ############
// Class header
// ############
class Rescale
{
 private:

  double q95_old; // Original q_95
  double q95_new; // New q_95
  double a1;      // Rescaling parameter

 public:

  // Constructor
  Rescale ();
  // Rescale equilibrium
  void RescaleEquilibrium ();
  // Rescale pFile
  void pFileRescale ();

  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open existing file for reading
  FILE* OpenFiler (char* filename);
  // Open existing file for appending
  FILE* OpenFilea (char* filename);
  // Call operating system 
  void CallSystem (char* command);  
};
  
#endif //RESCALE
