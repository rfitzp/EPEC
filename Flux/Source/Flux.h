// Flux.h

// ############################################################################
// Class to input gFile equilibrium data and output perturbed equilibrium data. 

// .................
// Calculation grid:
// .................

// Radial grid in PsiN = 1. - Psi/Psi_axis (assuming Psi = 0 on boundary) is 
//
// PsiN_j = PSILIM * (1. - (1. - s)^PACK)  for j = 0, NPSI-1
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

// All lengths normalized to R_0.
// All magnetic field-stengths normalized to B_0.

// ...................
// Inputs and outputs:
// ...................
// Calculation control parameters in namelist file Inputs/Flux.nml

// Equilibrium in Inputs/gFile or Inputs/gFiles
// Intermediate data in Outputs/Stage1.nc
// Final data in Outputs/Stage2.nc
// Data passed to programs NEOCLASSICAL and PHASE output to Outputs/fFile or Outputs/fFiles

// .........
// Versions:
// .........

// 1.0  - Initial version
// 1.2  - Improved gFile indexing
// 1.3  - Major rearrangement of input and output files
// 1.4  - Added linear interpolation
// 1.5  - Added RP1, Bt, Bt1, Bp, Bp1, and K_theta
// 1.6  - Removed QFLAG and Q95 functionality. Added calculation of A1, A2, A3 parameters
// 1.7  - Improved system calls
// 1.8  - Renamed Namelist. Removed A2 and A3 parameters (too much noise)
// 1.9  - Added PSILIM to fFile
// 1.10 - Added A2res, PSIPED, Pped to fFile
// 1.11 - Added -P option
// 1.12 - Included ZOFF parameter in gFiles
// 1.13 - Restrict program to only calculated q-profile for PsiN < PSIRAT
// 1.14 - Added q_hat calculation
// 1.15 - Added C1 and C2 calculation
// 1.16 - Added E, F, H calculation
// 1.17 - Upgraded to gsl_odeiv2. Better Cnc calculation.
// 1.18 - Added RMP coil calculation
// 1.19 - Removed RMP coil calculation and E vectors. 
//         Added NECTCDF output. Adapted to run with OMFIT.

// 2.0  - Completely switched to OMFIT mode
// 2.1  - Added cylindrical tearing mode calculation
// 2.2  - Removed command line options
// 2.3  - Added resistive wall
// 2.4  - Added cylindrical quantities to NETCDF file

// #####################################################################################

#ifndef FLUX
#define FLUX

#define MAXFILENAMELENGTH 500

#define VERSION_MAJOR 2
#define VERSION_MINOR 4

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include <blitz/array.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_gamma.h>
#include <netcdf.h>

using namespace blitz;

// Pointers to right-hand side function for adaptive integration
extern "C" int pRhs1 (double, const double[], double[], void*);
extern "C" int pRhs2 (double, const double[], double[], void*);
extern "C" int pRhs3 (double, const double[], double[], void*);
extern "C" int pRhs4 (double, const double[], double[], void*);
extern "C" int pRhs5 (double, const double[], double[], void*);
extern "C" int pRhs6 (double, const double[], double[], void*);
extern "C" int pRhs7 (double, const double[], double[], void*);
extern "C" int pRhs8 (double, const double[], double[], void*);

// Namelist reading function
extern "C" void NameListRead (int* INTG, int* NPSI, double* PACK, int* NTHETA, int* NNC, int* NTOR, double* H0,
			      double* ACC, double* ETA, int* MMIN, int* MMAX, double* PSILIM, double* TIME,
			      double* PSIPED, int* NSMOOTH, double* PSIRAT, int* NEOANG, double* RW);

// gFile reading function
extern "C" void gFileRead ();

// gFile interpolation functions
extern "C" void gFileInterpolateLinear    ();
extern "C" void gFileInterpolateQuadratic ();
extern "C" void gFileInterpolateCubic     ();
extern "C" void gFileInterpolateQuartic   ();

// ############
// Class header
// ############
class Flux
{
 private:

  // Control parameters read from Inputs/Flux.nml
  int    NTOR;    // Toroidal mode number
  int    MMIN;    // Minimum poloidal mode number
  int    MMAX;    // Maximum poloidal mode number

  double PSILIM;  // Maximum value of PsiN for safety-factor calculation
  double PSIRAT;  // All rational surfaces lying in region PsiN > PSIRAT are ignored
  double PSIPED;  // PsiN value at top of pedestal

  int    INTG;    // If != 0 then use interpolated gFile
  double TIME;    // Experimental time (ms)

  double RW;      // Radius of resistive wall (units of minor radius)

  int    NPSI;    // Number of points in PsiN grid
  double PACK;    // Packing index for PsiN grid
  int    NTHETA;  // Number of points in theta grid
  int    NNC;     // Number of neoclassical harmonics
  int    NSMOOTH; // Number of smoothing cycles for higher derivatives of q

  double H0;      // Initial integration step-length for equilibirum flux surface integrals 
  double ACC;     // Integration accuracy for equilibrium flux surface integrals
  double ETA;     // Regularization factor for Green's function
  int    NEOANG;  // Flag for using neoclassical angle in E-matrix calculation

  double EPS;     // Cylindrical tearing solutions launched from magnetic axis at r = EPS 
  double DELTA;   // Cylindrical tearing solutions integrated to r = r_s +/- DELTA 
  
  // Stage 1 parameters
  double          R0;       // Scale major radius (m)
  double          B0;       // Scale toroidal magnetic field strength (T)
  double          RLEFT;    // Bounding box coordinate
  double          ZLOW;     // Bounding box coordinate
  double          RRIGHT;   // Bounding box coordinate
  double          ZHIGH;    // Bounding box coordinate
  double          Raxis;    // Magnetic axis coordinate
  double          Zaxis;    // Magnetic axis coordinate
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
  double*         G;        // g(Psi)
  double*         Pr;       // p(Psi)
  double*         GGp;      // g dg/dPsi
  double*         Prp;      // dp/dPsi
  double*         Q;        // q(Psi)
  double*         Jphi;     // J_phi(Psi)

  double          dR2, dZ2, dR3, dZ3;
 
  // Flux coordinate construction parameters
  double  Psic;        // Psi on magnetic axis
  double  Rbound;      // R coordinate of plasma boundary on inboard midplane
  double  Rbound1;     // R coordinate of plasma boundary on outboard midplane
  int     ia;          // R grid index of inboard plasma boundary at Z = Z_axis
  int     ic;          // R grid index of magnetic axis
  int     jc;          // Z grid index of magnetic axis
  int     L;           // Number of points in Psi(R,Zaxis) array
  double* s;           // Array of s = sqrt[1 - Psi(R,Zaxis)] values
  double* Rs;          // Array of R(s) values, where R is major radius on inboard midplane
  double  q95;         // Safety-factor on 95% flux surface
  double  r95;         // Radial coordinate of 95% flux surface
  double  qrat;        // Safety-factor at PsiN = PSIRAT flux surface
  double  rrat;        // Radial coordinate at PsiN = PSIRAT flux surface
  double  qa;          // Safety-factor at plasma boundary
  double  ra;          // Radial coordinate of plasma boundary
  double  Pped;        // Pedestal pressure / central pressure

  double  qgp, qgp1, qgp2;
  int     ires;

  // Stage2 profile parameters
  double* P;           // Psi array
  double* RP;          // R(Psi) on inboard midplane
  double* RP1;         // R(Psi) on outboard midplane 
  double* Bt;          // B_toroidal(Psi) on inboard midplane
  double* Bt1;         // B_toroidal(Psi) on outboard midplane
  double* Bp;          // B_poloidal(Psi) on inboard midplane
  double* Bp1;         // B_poloidal(Psi) on outboard midplane
  double* rP;          // r(Psi)
  double* GP;          // g(Psi)
  double* QGP;         // q(Psi)/g(Psi) 
  double* QP;          // q(Psi)
  double* PP;          // P(Psi)
  double* GPP;         // dg/dPsi
  double* PPP;         // dP/dPsi
  double* JP;          // J(Psi)
  double* S;           // sqrt(1 - Psi)
  double* QX;          // q(Psi) from gFile
  double* J0;          // GGJ integral

  double* PsiN;        // PsiN array
  double* QPN;         // dQ/dPsiN array
  double* QPPN;        // d^2Q/dPsiN^2 array
  double* A1;          // QP/QPN/fabs(Psic) array
  double* A2;          // QPPN/QPN/3 array
  double* JPr;         // dJ/dr
  double* JPrr;        // d^2J/dr^2

  // Rational surface data
  int     nres;        // Number of rational surfaces
  int*    mres;        // Poloidal mode numbers at rational surfaces
  double* qres;        // Safety-factors at rational surfaces
  double* PsiNres;     // Normalized poloidal fluxes at rational surfaces   
  double* rres;        // Minor radii of rational surfaces
  double* sres;        // Magnetic shears at rational surfaces
  double* gres;        // g values at rational surfaces
  double* Rres;        // R coordinates of rational surfaces on inboard midplane
  double* Rres1;       // R coordinates of rational surfaces on outboard midplane
  double* gmres;       // gamma values on rational surfaces
  double* Ktres;       // K_t values at rational surfaces
  double* Kares;       // K_ast values at rational surfaces
  double* Khres;       // K_theta values at rational surfaces
  double* fcres;       // Fraction of circulating particles at rational surfaces
  double* ajj;         // Metric elements at rational surfaces
  double* dPsidr;      // dPsi/dr at rational surfaces
  double* q_hat;       // q_hat at rational surfaces
  double* A1res;       // A1 values at rational surfaces
  double* A2res;       // A2 values at rational surfaces
  double* C1res;       // C1 values at rational surfaces
  double* C2res;       // C2 values at rational surfaces
  double* alpres;      // alpha values at rational surfaces
  double* betres;      // beta values at rational surfaces
  double* gamres;      // gamma' values at rational surfaces

  // Glasser-Greene-Johnson data
  double* J1;          // GGJ integral
  double* J2;          // GGJ integral
  double* J3;          // GGJ integral
  double* J4;          // GGJ integral
  double* J5;          // GGJ integral
  double* J6;          // GGJ integral
  double* E;           // GGJ index
  double* F;           // GGJ index
  double* H;           // GGJ index
  
  // Straight angle flux coordinate data
  double*     th;      // theta array
  gsl_matrix* Rst;     // R versus theta on rational surfaces
  gsl_matrix* Zst;     // Z versus theta on rational surfaces

  // Neoclassical angle flux coordinate data
  double*     Th;      // Theta array
  gsl_matrix* theta;   // theta evaluated on Theta array on rational surfaces
  gsl_matrix* Rnc;     // R versus Theta on rational surfaces
  gsl_matrix* Znc;     // Z versus Theta on rational surfaces
  gsl_matrix* Bnc;     // B versus Theta on rational surfaces
  gsl_matrix* Cnc;     // dB/dTheta versus Theta on rational surfaces
  gsl_matrix* factor;  // Transformation function for neoclassical E-matrix interagration

  // Neoclassical parameter data
  double*     I1;      // Neoclassical integral
  double*     I2;      // Neoclassical integral
  double*     I3;      // Neoclassical integral
  gsl_matrix* I4;      // Neoclassical integrals
  gsl_matrix* I5;      // Neoclassical integrals
  gsl_matrix* I6;      // Neoclassical integrals
  double*     I7;      // Neoclassical integral
  double*     I8;      // Neoclassical integral

  // Tearing stability matrices
  gsl_matrix_complex* FF;  // F-matrix
  gsl_matrix_complex* EE;  // E-matrix

  // Cylindrical tearing stability data
  double* A;         // Island saturation parameters
  double* B;         // Island saturation parameters
  double* Delta_nw;  // Even parity tearing stability indices with no wall
  double* Sigma_nw;  // Odd parity tearing stability indices with no wall
  double* Delta_pw;  // Even parity tearing stability indices with perfect wall
  double* Sigma_pw;  // Odd parity tearing stability indices with perfect wall
  double* Delta_w;   // Wall stability index
  double* Sigma_w;   // Plasma/wall coupling index
  double* Poem1;     // Small-island POEM logarithmic saturation terms
  double* Poem2;     // Small-island POEM saturation terms
  double* Poem3;     // Large-island POEM saturation terms

  // Weights for Simpson's rule
  double          hh;
  Array<double,1> Weight1D;
  Array<double,1> weight1D;
  Array<double,2> Weight2D;
  
public:

  // Constructor
  Flux ();
  // Solve problem
  void Solve ();

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
  // Evaluate right-hand sides of GGJ equation
  int Rhs6 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of tearing equation inside rational surface
  int Rhs7 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of tearing equation outside rational surface
  int Rhs8 (double r, const double y[], double dydr[], void*);

private:

  // Set global parameters
  void SetParameters ();
  // Input gFile data and output Stage1 data
  void Stage1 ();
  // Input Stage1 data and output Stage2 data
  void Stage2 ();
  // Write Stage2 NETCDF file
  void WriteStage2Netcdfc ();

  // Set Simpson weights
  void Stage2SetSimpsonWeights ();
  // Read data for Stage2 calculations
  void Stage2ReadData ();
  // Calculate Stage2 q profile
  void Stage2CalcQ ();
  // Find Stage2 rational surfaces
  void Stage2FindRational ();
  // Calculate GGJ data at rational surfaces
  void Stage2CalcGGJ ();
  // Calculate Stage2 straight angle data at rational surfaces
  void Stage2CalcStraightAngle ();
  // Calculate Stage2 neoclassical angle data at rational surfaces
  void Stage2CalcNeoclassicalAngle ();
  // Calculate Stage2 neoclassical parameters at rational surfaces
  void Stage2CalcNeoclassicalPara ();
  // Calculate Stage2 tearing stability matrices
  void Stage2CalcMatrices ();
  // Calculate Stage2 cylindrical tearing mode indices
  void Stage2CalcTearing ();

  // Calculate q(P)/g(P) profile
  void CalcQGP ();
  // Check q value at rational surfaces
  void CheckQP ();
  // Calculate r(P) profile
  void CalcrP ();
  // Calculate GGJ parameters at rationa surface
  void CalcGGJ ();
  // Calculate straight angle data on rational surfaces
  void CalcStraightAngle ();
  // Calculate gamma values at rational surfaces
  void CalcGamma ();
  // Calculate neoclassical angle data on rational surfaces
  void CalcNeoclassicalAngle ();
  // Calculate no wall cylindrical tearing solutions
  void CalcTearingSolutionNoWall (int i);
  // Calculate perfect wall cylindrical tearing solutions
  void CalcTearingSolutionPerfectWall (int i);

  // gFile interpolation routines
  void gFileInterp          (vector<string> gFileName,   vector<double> gFileTime,   int gFileNumber, double time);
  void gFileInterpLinear    (char* gFile1, double time1, char* gFile,  double time);
  void gFileInterpQuadratic (char* gFile1, double time1, char* gFile2, double time2, char* gFile,     double time);
  void gFileInterpCubic     (char* gFile1, double time1, char* gFile2, double time2, char* gFile3,    double time3,
			     char* gFile,  double time);
  void gFileInterpQuartic   (char* gFile1, double time1, char* gFile2, double time2, char* gFile3,    double time3,
			     char* gFile4, double time4, char* gFile,  double time);
    
  // 1D interpolation function with nonuniform grid
  double Interpolate         (int I, double* X, double* Y, double x,                                 int order);
  double InterpolateCubic    (       double* X, double* Y, double x, int i0, int i1, int i2,         int order);
  double InterpolateQuartic  (       double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

  // 1D interpolation function with nonuniform grid for interpolating q/g profile
  double InterpolateQ          (int I, double* X, double* Y, double x,                                 int order, int method);
  double InterpolateQQuadratic (       double* X, double* Y, double x, int i0, int i1,                 int order);
  double InterpolateQCubic     (       double* X, double* Y, double x, int i0, int i1, int i2,         int order);
  double InterpolateQQuartic   (       double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

  // 1D interpolation function with nonuniform grid and periodic function
  double InterpolatePeriodic         (int I, double* X, double* Y, double x,                                 int order);
  double InterpolatePeriodicQuartic  (int I, double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

  // Interpolate Psi on uniform 2D grid
  double InterpolatePsi               (double RR, double ZZ,                                                                 int order);
  double InterpolatePsiCubicCubic     (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2,                 int order);
  double InterpolatePsiQuarticCubic   (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2,         int order);
  double InterpolatePsiCubicQuartic   (double RR, double ZZ, int i0, int i1, int i2, int j0, int j1, int j2, int j3,         int order);
  double InterpolatePsiQuarticQuartic (double RR, double ZZ, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order);

  // Evaluate dPsi/dR (R, Z)
  double GetPsiR (double r, double z);
  // Evaluate dPsi/dZ (R, Z)
  double GetPsiZ (double r, double z);
  // Evaluate d^2Psi/dR^2 (R, Z)
  double GetPsiRR (double r, double z);
  // Evaluate d^2Psi/dRdZ (R, Z)
  double GetPsiRZ (double r, double z);
  // Evaluate d^2Psi/dZ^2 (R, Z)
  double GetPsiZZ (double r, double z);

  // Savitsky-Gorlay smoothing algorithm
  void Smoothing (int N, double *y);
   
  // Calculate Green's functions
  double GreenPlasmaCos    (int i, int j, int ip, int jp);
  double GreenPlasmaSin    (int i, int j, int ip, int jp);
  double GreenPlasmaCosNeo (int i, int j, int ip, int jp);
  double GreenPlasmaSinNeo (int i, int j, int ip, int jp);

  // Calculate toroidal P function
  double ToroidalP (int m, int n, double z);

  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open existing file for reading
  FILE* OpenFiler (char* filename);
  // Open existing file for appending
  FILE* OpenFilea (char* filename);
  // Call operating system 
  void CallSystem (char* command);
};

#endif //FLUX
