// Flux.h

// ############################################################################
// Class to input gFile equilibrium data and output perturbed equilibrium data. 
//
// Radial grid in PsiN = 1. - Psi/Psi_axis (assuming Psi = 0 on boundary) is 
//
// PsiN_j = s  for j = 0, NPSI-1
//
//  where s = j /(NPSI-1),
//
// Poloidal grid in PEST poloidal angle theta is
//
// theta_k = 2*M_PI * t  for k = 0, NTHETA-1
//
//  where t = k /(NTHETA-1).
//
// theta = 0 on inboard midplane.
// theta > 0 above midplane.

// Command line options:
// -g INTG   - override INTG value from namlist
// -n NTOR   - override NTOR value from namelist
// -m MMIN   - override MMIN value from namelist
// -M MMAX   - override MMAX value from namelist
// -p PSILIM - override PSILIM value from namelist
// -t TIME   - sets experimental time

// Calculation control parameters in namelist file INPUTS/Flux.in.

// Equilibrium in Inputs/gFile.
// Intermediate data in folder Outputs/Stage1/
// Final data in folder Outputs/Stage2/
// Data passed to programs NEOCLASSICAL and PHASE output to Outputs/fFile

// Version:
// 1.0 - Initial version
// 1.2 - Improved gFile indexing
// 1.3 - Major rearrangement of input and output files

// ############################################################################

#ifndef FLUX
#define FLUX

#define VERSION_MAJOR 1
#define VERSION_MINOR 3

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <blitz/array.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_gamma.h>

#define MAXFILENAMELENGTH 500

using namespace blitz;

// Pointers to right-hand side function for adaptive integration
extern "C" int pRhs1 (double, const double[], double[], void*);
extern "C" int pRhs2 (double, const double[], double[], void*);
extern "C" int pRhs3 (double, const double[], double[], void*);
extern "C" int pRhs4 (double, const double[], double[], void*);
extern "C" int pRhs5 (double, const double[], double[], void*);

// Namelist reading function
extern "C" void NameListRead (int* INTG, int* NPSI, int* NTHETA, int* NNC, int* NTOR, int* QFLG, double* Q95, double* H0,
			      double* ACC, double* ETA, double* DR, int* MMIN, int* MMAX, double* PSILIM, double* TIME);

// gFile reading function
extern "C" void gFileRead ();

// gFile interpolation functions
extern "C" void gFileInterpolateQuadratic ();
extern "C" void gFileInterpolateCubic ();
extern "C" void gFileInterpolateQuartic ();

// ############
// Class header
// ############
class Flux
{
 private:
  
  // Control parameters read from Inputs/Flux.in
  int    NPSI;         // Number of points in PsiN grid
  int    NTHETA;       // Number of points in theta grid
  int    NNC;          // Number of neoclassical harmonics
  int    QFLG;         // QFLG = 0 - use q95 from gFile
                       // QFLG = 1 - rescale q95 to Q95 by adding constant to gg'
  double Q95;          // Target Q95 (QFLG = 1)
  int    NTOR;         // Toroidal mode number
  int    MMIN;         // Minimum poloidal mode number
  int    MMAX;         // Maximum poloidal mode number
  double PSILIM;       // Maximum PsiN for rational surface
  double H0;           // Initial integration step-length for equilibirum flux surface integrals 
  double ACC;          // Integration accuracy for equilibrium flux surface integrals
  double ETA;          // Regularization factor for Green's function
  double DR;           // Discritization parameter for simulated Mirnov data
  double TIME;         // Experimental time
  int    INTG;         // If != 0 then use interpolated gFile

  // Toroidal Mirnov coil array locations
  double RIN;          // Limiter R coordinate at inboard miplane
  double ROUT;         // Limiter R coordinate at outboard miplane
  
  // Stage 1 parameters
  double          R0;       // Scale major radius
  double          B0;       // Scale toroidal magnetic field strength
  int             NRPTS;    // Number of R points
  double*         RPTS;     // R array
  int             NZPTS;    // Number of Z points
  double*         ZPTS;     // Z array
  double          dR;       // R grid spacing
  double          dZ;       // Z grid spacing
  int             NBPTS;    // Number of boundary points
  double*         RBPTS;    // R on plasma boundary
  double*         ZBPTS;    // Z on plasma boundary
  int             NLPTS;    // Number of limiter points
  double*         RLPTS;    // R on limiter boundary
  double*         ZLPTS;    // Z on limiter boundary
  Array<double,2> PSIARRAY; // Psi (R, Z) array 
  double*         PSIN;     // PsiN array 
  double*         G ;       // g(Psi)
  double*         Pr;       // p(Psi)
  double*         GGp;      // g dg/dPsi
  double*         Prp;      // dp/dPsi
  double*         Q;        // q(Psi)

  double dR2, dZ2, dR3, dZ3;
 
  // Flux coordinate construction parameters
  double  Psic;        // Psi on magnetic axis
  double  Raxis;       // R coordinate of magnetic axis
  double  Zaxis;       // Z coordinate of magnetic axis
  double  Rbound;      // R coordinate of plasma boundary on inboard midplane
  int     ia;          // R grid index of inboard plasma boundary at Z=Z_axis
  int     ic;          // R grid index of magnetic axis
  int     jc;          // Z grid index of magnetic axis
  int     L;           // Number of points in Psi(R,Zaxis) array
  double* s;           // Array of s = sqrt[1 - Psi(R,Zaxis)] values
  double* Rs;          // Array of R(s) values
  double  q95;         // Safety-factor on 95% flux surface
  double  r95;         // Radial coordinate of 95% flux surface
  double  qlim;        // Safety-factor at PsiN = PSILIM flux surface
  double  rlim;        // Radial coordinate at PsiN = PSILIM flux surface
  double  qa;          // Safety-factor at plasma boundary
  double  ra;          // Radial coordinate of plasma boundary
  double  qgp, qgp1;

  // Stage2 profile parameters
  double* P ;          // Psi array
  double* RP;          // R(Psi)
  double* rP;          // r(Psi)
  double* GP;          // g(Psi)
  double* QGP;         // q(Psi)/g(Psi) 
  double* QP;          // q(Psi)
  double* PP;          // P(Psi)
  double* GPP;         // dg/dPsi
  double* PPP;         // dP/dPsi
  double* S;           // sqrt(1 - Psi)
  double* QX;          // q(Psi) from gFile

  // Rational surface data
  int     nres;        // Number of rational surfaces
  int*    mres;        // Poloidal mode numbers at rational surfaces
  double* qres;        // Safety-factors at rational surfaces
  double* PsiNres;     // Normalized poloidal fluxes at rational surfaces   
  double* rres;        // Minor radii of rational surfaces
  double* sres;        // Magnetic shears at rational surfaces
  double* gres;        // g values at rational surfaces
  double* Rres;        // R coordinates of rational surfaces on midplane
  double* gmres;       // gamma values on rational surfaces
  double* Ktres;       // K_t values at rational surfaces
  double* Kares;       // K_ast values at rational surfaces
  double* fcres;       // Fraction of circulating particles at rational surfaces
  double* ajj;         // Metric elements at rational surfaces
  double* dPsidr;      // dPsi/dr at rational surfaces
  
  // Straight angle flux coordinate data
  double*     th;      // theta array
  gsl_matrix* Rst;     // R versus theta on rational surfaces
  gsl_matrix* Zst;     // Z versus theta on rational surfaces

  // Neoclassical angle flux coordinate data
  double*     Th;      // Theta array
  gsl_matrix* Rnc;     // R versus Theta on rational surfaces
  gsl_matrix* Znc;     // Z versus Theta on rational surfaces
  gsl_matrix* Bnc;     // B versus Theta on rational surfaces
  gsl_matrix* Cnc;     // dB/dTheta versus Theta on rational surfaces

  // Neoclassical parameter data
  double*     I1;      // Neoclassical integral
  double*     I2;      // Neoclassical integral
  double*     I3;      // Neoclassical integral
  gsl_matrix* I4;      // Neoclassical integrals
  gsl_matrix* I5;      // Neoclassical integrals
  gsl_matrix* I6;      // Neoclassical integrals

  // Perturbed equilibrium data
  gsl_matrix_complex* FF;  // F-matrix
  gsl_matrix_complex* EE;  // E-matrix
  gsl_vector_complex* EI;  // Response vector for inboard toroidal Mirnov array
  gsl_vector_complex* EO;  // Response vector for outboard toroidal Mirnov array

  // Weights for Simpson's rule
  double hh;
  Array<double,1> Weight1D;
  Array<double,1> weight1D;
  Array<double,2> Weight2D;
  
public:

  // Constructor
  Flux ();
  // Solve problem
  void Solve (int _INTG, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM);

  // Evaluate right-hand sides of q/g equation
  int Rhs1 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of r equation
  int Rhs2 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of theta equation
  int Rhs3 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of gamma equation
  int Rhs4 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of Theta equation
  int Rhs5 (double r, const double y[], double dydr[], void*);

private:

  // Set global parameters
  void SetParameters (int _INTG, int _NTOR, int _MMIN, int _MMAX, double _TIME, double _PSILIM);
  // Input gFile data and output Stage1 data
  void Stage1 ();
  // Input Stage1 data and output Stage2 data
  void Stage2 ();

  // Set Simpson weights
  void Stage2SetSimpsonWeights ();
  // Read data for Stage2 calculations
  void Stage2ReadData ();
  // Calculate Stage2 q profile
  void Stage2CalcQ ();
  // Find Stage2 rational surfaces
  void Stage2FindRational ();
  // Calculate Stage2 straight angle data at rational surfaces
  void Stage2CalcStraightAngle ();
  // Calculate Stage2 neoclassical angle data at rational surfaces
  void Stage2CalcNeoclassicalAngle ();
  // Calculate Stage2 neoclassical parameters at rational surfaces
  void Stage2CalcNeoclassicalPara ();
  // Calculate Stage2 stability matrices
  void Stage2CalcMatrices ();

  // Calculate q(P)/g(P) profile
  void CalcQGP ();
  // Check q value at rational surfaces
  void CheckQP ();
  // Calculate r(P) profile
  void CalcrP ();
  // Calculate straight angle data on rational surfaces
  void CalcStraightAngle ();
  // Calculate gamma values at rational surfaces
  void CalcGamma ();
  // Calculate neoclassical angle data on rational surfaces
  void CalcNeoclassicalAngle ();

  // gFile interpolation routines
  void gFileInterp          (vector<string> gFileName,   vector<double> gFileTime,   int gFileNumber, double time);
  void gFileInterpQuadratic (char* gFile1, double time1, char* gFile2, double time2, char* gFile,     double time);
  void gFileInterpCubic     (char* gFile1, double time1, char* gFile2, double time2, char* gFile3,    double time3,
			     char* gFile,  double time);
  void gFileInterpQuartic   (char* gFile1, double time1, char* gFile2, double time2, char* gFile3,    double time3,
			     char* gFile4, double time4, char* gFile,  double time);
    
  // 1D interpolation function with nonuniform grid
  double Interpolate         (int I, double* X, double* Y, double x,                                 int order);
  double InterpolateCubic    (       double* X, double* Y, double x, int i0, int i1, int i2,         int order);
  double InterpolateQuartic  (       double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

  // 1D interpolation function with nonuniform grid and periodic function
  double InterpolatePeriodic         (int I, double* X, double* Y, double x,                                 int order);
  double InterpolatePeriodicQuartic  (int I, double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

  // Interpolate Psi on uniform 2D grid
  double InterpolatePsi               (double RR, double ZZ,                                                                 int order);
  double InterpolatePsiCubicCubic     (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2,                 int order);
  double InterpolatePsiQuarticCubic   (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2,         int order);
  double InterpolatePsiCubicQuartic   (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2, int j3,         int order);
  double InterpolatePsiQuarticQuartic (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order);

  // Evaluate Psi (R, Z)
  double GetPsi (double r, double z);
  // Evaluate dPsi/dR (R, Z)
  double GetPsiR (double r, double z);
  // Evaluate dPsi/dZ (R, Z)
  double GetPsiZ (double r, double z);
 
  // Calculate Green's functions
  double GreenPlasmaCos   (int i, int j, int ip, int jp);
  double GreenPlasmaSin   (int i, int j, int ip, int jp);
  double GreenInboardCos  (int i, int j);
  double GreenInboardSin  (int i, int j);
  double GreenOutboardCos (int i, int j);
  double GreenOutboardSin (int i, int j);

  // Calculate toroidal P function
  double ToroidalP (int m, int n, double z);

  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open existing file for reading
  FILE* OpenFiler (char* filename);
  // Open existing file for appending
  FILE* OpenFilea (char* filename);
};

#endif //FLUX
