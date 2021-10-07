// Rescale.h

// ##############################################################
// Class to rescale equilibrium gFile, pFile, and cFile

// Type 1:
// Rescaling of density by factor SCALE (from namelist)

// Type 2:
// Rescaling of temperature by factor SCALE (from namelist)

// Type 3:
// Rescaling of plasma size by factor SCALE (from namelist)

// Type 4:
// Rescaling of diffusivities by factor SCALE (from namelist)

// Type 5:
// Shift of pressure profile by PSHIFT (kPa) (from namelist)

// Type 6:
// Shift of ExB frequency profile by WSHIFT (krad/s) (from namelist)

// Type 7:
// Shift of safety factor profile such that q(PsiN=0.95) = Q_95 (from namelist)

// Initial gFile in Inputs/gFile
// Rescaled gFile in Outputs/gFile

// Initial pFile in Inputs/pFile
// Rescaled pFile in Outputs/pFile

// Initial cFile in Inputs/cFile
// Rescaled cFile in Outputs/cFile

// 1.0 - Original version
// 2.0 - Extended version

// ##############################################################

#ifndef RESCALE
#define RESCALE

#define VERSION_MAJOR 2
#define VERSION_MINOR 0
#define MAXFILENAMELENGTH 500

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include "Field.h"

// Read namelist
extern "C" void NameListRead (int* TYPE, double* SCALE, double* PSHIFT, double* WSHIFT, double* Q95, int* OPOINT, int* XPOINT);

// Perform Type 0 gFile rescaling 
extern "C" void gFileRescaleType0 ();

// Perform Type I gFile rescaling 
extern "C" void gFileRescaleTypeI (double* A, int* OPOINT, int* XPOINT);

// Perform Type II gFile rescaling
extern "C" void gFileRescaleTypeII (double* A, int* OPOINT, int* XPOINT);

// Perform Type III gFile rescaling 
extern "C" void gFileRescaleTypeIII (double* PSHIFT);

// Perform Type 7 gFile rescaling 
extern "C" void gFileRescaleType7 (double* Q95, int* OPOINT, int* XPOINT, double* q95_old, double* a1);

// ############
// Class header
// ############
class Rescale
{
 private:

  int    TYPE;    // Rescaling type (1=density, 2=temperature, 3=size, 4=diffusivity, 5=pressure, 6=ExB freqeuncy, 7=safety factor)
  double SCALE;   // Rescaling factor (rescaling types 1-4)
  double PSHIFT;  // Pressure profile shift factor (kPa) (rescaling type 5)
  double WSHIFT;  // ExB frequency profile shift factor (krad/s) (rescaling type 6)
  double Q95;     // Target safety-factor at 95% flux-surface (rescaling type 7)
  int    OPOINT;  // Flag for finding new O-point
  int    XPOINT;  // Flag for finding new X-point

  double a1;      // Pressure rescaling parameter (rescaling type 7)
  
 public:

  // Constructor
  Rescale ();

  // Rescale plasma equilibrium
  void RescaleEquilibrium ();

  // Extract 1/2/ne and 1/2/ni fields from pFile
  void GetnFields (Field& ine, Field& ini);
  // Extract RBp and R fields from pFile
  void GetRFields (Field& RBp, Field& R);

  // Perform Type 1 rescaling of pFile
  void pFileRescaleType1 (double an); 
  // Perform Type 2 rescaling of pFile
  void pFileRescaleType2 (double at);
  // Perform Type 3 rescaling of pFile
  void pFileRescaleType3 (double ar);
  // Perform Type 4 rescaling of pFile
  void pFileRescaleType4 (double ac);
  // Perform Type 5 rescaling of pFile 
  void pFileRescaleType5 (double aw, Field& ine, Field& ini);
  // Perform Type 6 rescaling of pFile 
  void pFileRescaleType6 (double aw, Field& RBp, Field& R);
  // Perform Type 7 rescaling of pFile
  void pFileRescaleType7 (double a1);

  // Perform Type 4 rescaling of cFile
  void cFileRescaleType4 (double ac);

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
