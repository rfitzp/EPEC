// Phase.h

// #######################################################################
// Class to calculate time evolution of magnetic island chains interacting
// with resonant magnetic perturbation in toroidal tokamak plasma.
// Reads data from programs FLUX, NEOCLASSICAL, and GPEC.

// Inputs and outputs:
// ...................
// Calculation control parameters in namelist file Inputs/Phase.nml
// RMP wavform data in Inputs/Waveform.nml

// FLUX data in Inputs/fFile or Inputs/fFiles
// NEOCLASSICAL data in Inputs/nFile or Inputs/nFiles
// GPEC data in Inputs/uFile, Inputs/mFile, Inputs/lFile or Inputs/uFiles, Inputs/mFiles, Inputs/lFiles

// Intermediate data in Stage4.nc
// Final data in Stage5.nc or Stage5/*.nc

// ............................
// Major stages in calculation:
// ............................

// Stage 4 - Class reads FLUX/NEOCLASSICAL/GPEC data and plots error-field
//           drive versus relative phases of RMP coil currents for 1kA currents

// Stage 5 - Class performs island dynamics simulation

// .........
// Versions:
// .........

// 1.0  - Initial version
// 1.1  - Improved indexing of fFiles, nFiles, uFiles, and lFiles
// 1.2  - Major rearrangement of input and output files
// 1.3  - Added PsiN and island width in PsiN to Stage6 output files
// 1.4  - Added linear iterpolation
// 1.5  - Added wnl
// 1.6  - Redefined Sk. Corrected composite factor in Rutherford equation.
// 1.7  - Added SCALE as input value

// 2.0  - Relaxed no-slip constraint
// 2.1  - Redefined Sk. Corrected composite factor in Rutherford equation.
// 2.2  - Added SCALE as input value
// 2.3  - Separated waveform input data from main input data.
//        Modified finite island-width natural frequency interpolation.
// 2.4  - Added LIN flag
// 2.5  - Added middle coil set
// 2.6  - Limited island width to stop them extending beyond neighbouring rational surfaces
// 2.7  - Improved calculation of island widths
// 2.8  - Replaced TIME by TSTART and TEND. Time now input in ms.
// 2.9  - Renamed Namelist. Modified island width calculation.
// 2.10 - Added CHIR parameter
// 2.11 - Added FREQ == 2 option
// 2.12 - Added total pressure decrement
// 2.13 - Normalized total pressure decrement by P(0). Added IRMP.
// 2.14 - Added higher order transport calculation
// 2.15 - Bug fixes. Improved uFile/mFile/lFile interpolation.
// 2.16 - Added COPT flag
// 2.17 - Changed definition of MID
// 2.18 - Added more COPT options
// 2.19 - Natural frequency limited to linear/nonlinear model
// 2.20 - Added FREQ flag
// 2.21 - Use q_hat instead of q in velocity evolution equation
// 2.22 - Added velocity changes
// 2.23 - Added charge exchange damping time into angular equations of motion
// 2.24 - Added bootstrap/curvature terms into Rutherford equations
// 2.25 - Added polarization terms into Rutherford equations
// 2.26 - Allowed for EPEC RMP coil data
// 2.27 - Removed EPEC RMP coil data. Added Netcdf output. Adapted for use with OMFIT.

// 3.0  - Completely switched to OMFIT mode
// 3.1  - Added POEM terms
// 3.2  - Removed command line options
// 3.3  - Added resistive wall
// 3.4  - Added TOFF and electromagnetic torques
// 3.5  - Added spike and repeated ramp waveforms

// #######################################################################

#ifndef PHASE
#define PHASE

#define VERSION_MAJOR 3
#define VERSION_MINOR 5
#define MAXFILENAMELENGTH 500
#define MAXCONTROLPOINTNUMBER 500
#define MAXULFILELINELENGTH 500

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <blitz/array.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <netcdf.h>

using namespace blitz;

// Namelist function
extern "C" void NameListRead (int* NFLOW, int* STAGE5, int* INTF, int* INTN, int* INTU, int* NATS, int* OLD, int* FREQ, int* LIN, int* MID, int* COPT,
			      double* DT, double* TSTART, double* TEND, double* TOFF, double* SCALE, double* PMAX, double* CHIR, int* HIGH, int* RATS,
			      double* CORE, double* FFAC, int* CXD, int* BOOT, int* CURV, int* POLZ, double* TAUW, 
			      int* TYPE, int* NCTRL, double* TCTRL, double* ICTRL, double* PCTRL,
			      double* SSTART, double* SEND, double* SAMP, double* SPHA, double* SPVE, double* BACK,
			      double* RPERIOD, double* RSTART, double* REND, double* RPHA);

// ############
// Class header
// ############
class Phase
{
 private:

  // ------------
  // Program data
  // ------------

  // Read from Inputs/Phase.nml
  double  PMAX;    // Stage 4 phase scan from - PMAX*M_PI/2 to +PMAX*M_PI/2
  int     NPHA;    // Number of points in phase scan

  int     STAGE5;  // If != 0 then Stage5 calculation performed, otherwise calculation terminates after Stage4
  double  TSTART;  // Simulation start time (ms)
  double  TEND;    // Simulation end time (ms)
  double  DT;      // Data recorded every DT seconds
  double  TOFF;    // If OLD == 0 then run simulation for TOFF ms before recording data

  int     MID;     // Number of RMP coil sets (1, 2, or 3)
                   //  If == 1 requires lFile
                   //  If == 2 requires lFile and uFile
                   //  If == 3 requires lFile, uFile, and mFile
  int     COPT;    // If == 0 then no coil current optimization
                   // If == 1 then coil currents optimized in restricted fashion to maximize drive at closest rational surface to pedestal top
                   // If == 2 then coil currents optimized in unrestricted fashion to maximize drive at closest rational surface to pedestal top
                   // If == 3 then coil currents optimized in unrestricted fashion to maximize drive at pedestal top and minimize drive in core
  double  CORE;    // Core drive minimization factor for COPT = 3 (0.0 = no minimization, 1.0 = complete minmization)
  int     FREQ;    // Natural frequency switch:
                   //  If == 0 then use linear/nonlinear natural frequency with linear layer width as switch
                   //  If == 1 then use linear/nonlinear natural frequency with pressure flattening width as switch 
                   //  If == 2 then use linear/ExB/nonlinear natural frequency 
                   //  If == 3 then w_natural = FFAC * w_linear + (1 - FFAC) * w_EB
  double  FFAC;    // Natural frequency parameter (for FREQ = 3)

  int     LIN;     // If != 0 then perform purely linear calculation
  int     CXD;     // If != 0 include charge exchange damping in plasma angular equations of motion
  int     BOOT;    // If != 0 include effect of bootstrap current in Rutherford equations
  int     CURV;    // If != 0 include effect of magnetic field-line curvature in Rutherford equations
  int     POLZ;    // If != 0 include effect of ion polarization current in Rutherford equations

  double  TAUW;    // Wall time constant (ms)
  
  double  CHIR;    // Maximum allowable Chirikov parameter for vacuum islands

  int     INTF;    // If != 0 then use interpolated fFile 
  int     INTN;    // If != 0 then use interpolated nFile
  int     INTU;    // If != 0 then use interpolated lFile, uFile, and mFile
  int     NATS;    // If != 0 then use linear-only nFile interpolation
  int     RATS;    // If != 0 use linear-only interpolation for uFiles/mFiles/lFiles
  int     OLD;     // If != 0 then restart old calculation

  int     HIGH;    // If != 0 use higher-order transport analysis
  double  SCALE;   // GPEC scalefactor
  int     NFLOW;   // Number of flow harmonics included in model

  int     TYPE;    // RMP waveform type (1=programmed, 2=spike, 3=repeated ramp) 

  int     NCTRL;   // Number of programmed RMP control points
  double* TCTRL;   // Programmed RMP control times (ms)
  double* ICTRL;   // Programmed peak current flowing in RMP coils at control times (kA) 
  double* PCTRL;   // Programmed relative phases of RMP coil currents (units of PI) at control times

  double  SSTART;  // Time of start of RMP spike (ms)
  double  SEND;    // Time of end of RMP spike (ms)
  double  WAMOD;   // Amplitude modulation frequency of RMP coil currents duing spike (krad/s)
  double  WPMOD;   // Realtive phase velocity of RMP coil currents during spike (krad/s)
  double  SAMP;    // RMP coils current amplitudes at start of spike (kA)
  double  SPHA;    // Relative phase of RMP coil currents at start of spike (units of PI)
  double  SPVE;    // Phase velocity of RMP coil currents during spike (krad/s)
  double  BACK;    // Background RMP coil currents (kA)

  double  RPERIOD; // Repeated RMP ramp repetition period (ms) 
  double  RSTART;  // RMP coil currents at start of ramp (first ramp starts at t=0) (kA)
  double  REND;    // RMP coil currents at end of ramp (kA)
  double  RPHA;    // Relative phase of RMP coil currents during ramp (units of PI)

  // ----------------------
  // Data from program FLUX
  // ----------------------

  // Read from Inputs/fFile
  double R_0;             // Scale major radius (m)
  double B_0;             // Scale toroidal field-strength (T)
  double q95;             // Safety factor at PsiN = 0.95 flux surface
  double r95;             // Normalized radius of PsiN = 0.95 flux surface
  double qrat;            // Safety factor at PsiN = PSIRAT flux surface
  double rrat;            // Normalized radius of PsiN = PSIRAT flux surface
  double q0;              // Safety factor at magnetic axis
  double qa;              // Safety factor at plasma boundary
  double PSILIM;          // Limiting value of PsiN
  double PSIPED;          // Value of PsiN at top of pedestal
  double PSIRAT;          // Value of PsiN beyond which resonant surfaces ignored
  double Pped;            // Pedestal pressure / central pressure
  int    nres;            // Number of resonant surfaces in plasma
  gsl_matrix_complex* FF; // Plasma inverse tearing stability matrix
  gsl_matrix_complex* EE; // Plasma tearing stability matrix 
 
  Array<double,2> FFh;    // Moduli of FF elements
  Array<double,2> EEh;    // Moduli of EE elements    
  Array<double,2> xih;    // Arguments of EE elements
  Array<double,1> A1;     // A1 values at rational surfaces
  Array<double,1> Poem1;  // POEM island saturation terms at rational surfaces
  Array<double,1> Poem2;  // POEM island saturation terms at rational surfaces
  Array<double,1> Poem3;  // POEM island saturation terms at rational surfaces
  Array<double,1> Deltaw; // Wall stability index
  Array<double,1> Sigmaw; // Plasma/wall coupling constant

  // ------------------------------
  // Data from program NEOCLASSICAL
  // ------------------------------
 
  // Read from Inputs/nFile
  double          tau_A;    // Alfven time (s)
  double          P0;       // Central thermal pressure (10^19 m^-3 keV)
  Array<int,1>    mk;       // Resonant poloidal mode numbers
  Array<int,1>    ntor;     // Resonant toroidal mode number
  Array<double,1> rk;       // Minor radii of resonant surfaces / a
  Array<double,1> qk;       // Safety-factors at resonant surfaces
  Array<double,1> rhok;     // Normalized mass densities at resonant surfaces
  Array<double,1> a;        // Normalized (to R_0) plasma minor radius
  Array<double,1> Sk;       // Lundquist numbers at resonant surfaces
  Array<double,1> chipk;    // Momentum diffusivities at resonant surfaces (m^2/s)
  Array<double,1> taumk;    // Normalized momentum confinement timescales at resonant surfaces
  Array<double,1> tautk;    // Normalized poloidal flow damping timescales at resonant surfaces
  Array<double,1> tauxk;    // Normalized charge exchange damping timesscales at resonant surfaces
  Array<double,1> fack;     // Island width factors at resonant surfaces
  Array<double,1> delk;     // Linear layer widths at resonant surfaces
  Array<double,1> wkl;      // Normalized linear natural frequencies at resonant surfaces
  Array<double,1> wke;      // Normalized ExB natural frequencies at resonant surfaces
  Array<double,1> wkn;      // Normalized nonlinear natural frequencies at resonant surfaces
  Array<double,1> etae;     // Relative electron temperature gradients at resonant surfaces
  Array<double,1> etai;     // Relative ion temperature gradients at resonant surfaces
  Array<double,1> dnedrk;   // Electron density gradients at resonant surfaces (in r) (10^19/m^-4)
  Array<double,1> dTedrk;   // Electron temperature gradients at resonant surfaces (in r) (keV/m)
  Array<double,1> Wcrnek;   // Critical island widths for density flattening at resonant surfaces (in r) (m)
  Array<double,1> WcrTek;   // Critical island widths for electron temperature flattening at resonant surfaces (in r) (m)
  Array<double,1> WcrTik;   // Critical island widths for ion temperature flattening at resonant surfaces (in r) (m)
  Array<double,1> akk;      // Metric elements at resonant surfaces
  Array<double,1> gk;       // g values at resonant surfaces
  Array<double,1> qhatk;    // q_hat values at resonant surfaces
  Array<double,1> C1k;      // C1 values at resonant surfaces
  Array<double,1> C2k;      // C2 values at resonant surfaces
  Array<double,1> alphabe;  // Electron bootstrap terms at resonant surfaces
  Array<double,1> alphabi;  // Ion bootstrap terms at resonant surfaces
  Array<double,1> alphac;   // Curvature terms at resonant surfaces
  Array<double,1> alphap;   // Ion polarization terms at resonant surfaces
  Array<double,1> rhothe;   // Electron poloidal gyroradii at resonant surfaces
  Array<double,1> rhothi;   // Ion poloidal gyroradii at resonant surfaces
  Array<double,1> dPsiNdr;  // R_0 dPsiN/dr values at resonant surfaces
  Array<double,1> PsiN;     // PsiN values at resonant surfaces
  Array<double,1> nek;      // Electron number densities at resonant surfaces (10^19/m^-3)
  Array<double,1> nik;      // Ion number densities at resonant surfaces (10^19/m^-3)
  Array<double,1> Tek;      // Electron temperatures at resonant surfaces (keV)
  Array<double,1> Tik;      // Ion temperatures at resonant surfaces (keV)
  Array<double,1> dnidrk;   // Ion density gradients at resonant surfaces (in r) (10^19/m^-4)
  Array<double,1> dTidrk;   // Ion temperature gradients at resonant surfaces (in r) (keV/m)
  Array<double,1> Factor1;  // Transport factor
  Array<double,1> Factor2;  // Transport factor
  Array<double,1> Factor3;  // Transport factor
  Array<double,1> Factor4;  // Transport factor
  Array<double,1> Factor5;  // Transport factor
  Array<double,1> Factor6;  // Transport factor
  Array<double,1> Factor7;  // Transport factor
  Array<double,1> Factor8;  // Transport factor
  Array<double,1> Factor9;  // Transport factor
  Array<double,1> Factor10; // Transport factor
  Array<double,1> Factor11; // Transport factor
  Array<double,1> Factor12; // Transport factor
  Array<double,1> Deltakp;  // Delta_k+ values at rational surfaces
  Array<double,1> Deltakm;  // Delta_k- values at rational surfaces
   
  // ----------------------
  // Data from program GPEC
  // ----------------------
  
  // Read from Inputs/uFile
  gsl_vector_complex* DeltaU; // Delta values for 1kA in upper RMP coil
  gsl_vector_complex* ChiU;   // Chi values for 1kA in upper RMP coil

  // Read from Inputs/mFile
  gsl_vector_complex* DeltaM; // Delta values for 1kA in middle RMP coil
  gsl_vector_complex* ChiM;   // Chi values for 1kA in middle RMP coil

  // Read from Inputs/lFile
  gsl_vector_complex* DeltaL; // Delta values for 1kA in lower RMP coil
  gsl_vector_complex* ChiL;   // Chi values for 1kA in lower RMP coil

  // ------------------
  // Stage 4 parameters
  // ------------------
  Array<double,1> phase;      // RMP coil phases  
  Array<double,2> wvac;       // Vacuum island widths (in PsiN) versus RMP coil phases

  // ----------------
  // Velocity factors
  // ----------------
  Array<double,1> j0p;   // Zeros of J0 Bessel function
  Array<double,1> j1p;   // Zeros of J1 Bessel function
  Array<double,2> torp;  // Velocity factors in poloidal equations of motion
  Array<double,2> tort;  // Velocity factors in toroidal equations of motion
  Array<double,3> natp;  // Poloidal velocity factors in natural frequency equations
  Array<double,3> natt;  // Toroidal velocity factors in natural frequency equations

  // ---------------
  // Simulation data
  // ---------------
  double          TIME;    // Time in ms
  double          Tstart;  // Normalized simulation start time
  double          Tend;    // Normalized simulation end time
  double          Toff;    // Normalized offset time
  Array<double,1> TT;      // Normalized type1 waveform control times
  double          dTT;     // Normalized recording time interval
  double          tauw;    // Normalized wall time constant
  double          Sstart;  // Normaized type 2 RMP spike start time
  double          Send;    // Normaized type 2 RMP spike end time
  double          Spve;    // Normalized type 2 RMP phase velocity
  double          Rperiod; // Normaized type 3 RMP repeated ramp repeat time
  double          irmp;    // Peak current flowing in RMP coils (kA) at current time
  double          prmp;    // Relative phases of RMP coil currents (radians) at current time
  Array<double,1> Psik;    // Normalized magnitutudes of reconnected fluxes at resonant surfaces
  Array<double,1> phik;    // Helical phases of reconnected fluxes at resonant surfaces
  Array<double,1> Xk;      // Cartesian components of reconnected fluxes at resonant surfaces
  Array<double,1> Yk;      // Cartesian components of reconnected fluxes at resonant surfaces
  Array<double,1> Psiwk;   // Normalized magnitutudes of fluxes at resonant wall
  Array<double,1> phiwk;   // Helical phases of fluxes at wall
  Array<double,1> Uk;      // Cartesian components of fluxes at wall
  Array<double,1> Vk;      // Cartesian components of fluxes at wall
  Array<double,2> alphakp; // Poloidal velocity factors
  Array<double,2> betakp;  // Toroidal velocity factors
  Array<double,1> ww;      // Island frequencies
  Array<double,1> vp;      // Island phase velocities
  Array<int,1>    lock;    // Locking flags
  Array<double,1> chi;     // RMP amplitudes at rational surfaces
  Array<double,1> zeta;    // RMP phases at rational surfaces

  // ------------------
  // Physical constants
  // ------------------
  double e;         // Magnitude of electron charge (SI)
  double epsilon_0; // Electric permittivity of free space (SI)
  double mu_0;      // Magnetic permeability of free space (SI)
  double m_p;       // Mass of proton (SI)
  double m_e;       // Mass of electron (SI)
  
  // --------------------------------
  // Adapative integration parameters
  // --------------------------------
  double acc;      // Integration accuracy
  double h0;       // Initial step-length
  double hmin;     // Minimum step-length
  double hmax;     // Maximum step-length
  int    maxrept;  // Maximum number of step recalculations
  double omegamax; // Maximum natural frequeny magnitude (krad/s)

  // ----
  // Misc
  // ----
  int count;
  
  // ----------------------
  // Public class functions
  // ----------------------
 public:

  // Constructor
  Phase ();
  // Destructor
  virtual ~Phase () {};  

  // Solve problem
  void Solve ();
  // Read data
  void Read_Data ();
  // Calculate vacuum flux versus relative phases of RMP coil currents
  void Scan_Shift ();
  // Calculate velocity factors
  void Calc_Velocity ();
  // Initialize island dynamics calculation
  void Initialize ();
  // Perform island dynamics calculation
  void IslandDynamics ();
  // Save island dynamics calculation
  void Save ();

  // Write Stage4 NETCDF file
  void WriteStage4Netcdfc ();

  // Calculate Irmp and Prmp
  void CalcRMP (double t);
  // Calculate coil currents and phases
  void CalcCoil (double t, double& IU, double& IM, double& IL, double& PU, double& PM, double& PL);
  // Find resonant surfaces that straddle top of pedestal
  int Findk ();
  // Find maximum of restricted three-coil function
  double FindMax (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML);
  // Find minimum of restricted three coil function
  double FindMin (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML);
  // Evaluate restricted three-coil function and its derivative
  void ThreeCoil (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML, double Delta, double& fun, double& deriv, double& dderiv);
  // Calculate chi and zeta vectors
  void CalcChiZeta (double t);
  // Pack simulation variables into single vector
  void Pack (Array<double,1> y);
  // Unpack simulation variables from single vector
  void Unpack (Array<double,1> y);
  // Pack right-hand sides into single vector
  void PackRhs (Array<double,1> XkRHS,      Array<double,1> YkRHS,
		Array<double,1> UkRHS,      Array<double,1> VkRHS,
		Array<double,2> alphakpRHS, Array<double,2> betakpRHS,
		Array<double,1> dydt);
  // Calculate natural frequency
  double GetNaturalFrequency (int j);
  // Calculate actual frequency
  double GetActualFrequency (int j);
  // Calculate change in natural frequency (krad/s)
  double GetDeltaOmega (int j);
  // Calculate change in poloidal angular velocity (krad/s)
  double GetDeltaOmegaTheta (int j);
  // Calculate change in toroidal angular velocity (krad/s)
  double GetDeltaOmegaPhi (int j);
  // Function to calculate change in toroidal velocity (km/s)
  double GetDeltaVPhi (int j);
  // Calculate change in parallel velocity (km/s)
  double GetDeltaVParallel (int j);
  // Calculate change in ExB velocity (km/s)
  double GetDeltaVEB (int j);
  // Calculate change in Er (kV/m)
  double GetDeltaEr (int j);
  // Calculate electromagnetic torques at rational surfaces
  void GetElectromagneticTorques (Array<double,1 >& y, double* T_Rmp, double* T_Wall, double* T_Tear);
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, Array<double,1>& y, Array<double,1>& dydx);
  // Adaptive-step integration routine
  void RK4Adaptive (double& x, Array<double,1>& y, double& h, double& t_err, 
		    double acc, double S, int& rept, int maxrept, 
		    double h_min, double h_max, int flag, int diag, FILE* file);
  // Fixed step integration routine
  void RK4Fixed (double& x, Array<double,1>& y, double h);

  // Interpolate uFiles
  void uFileInterp               (vector<string> uFileName,   vector<double> uFileTime,   int uFilenumber, double TIME);
  void uFileInterpolateLinear    (char* uFile1, double time1, char* uFile,  double time);
  void uFileInterpolateQuadratic (char* uFile1, double time1, char* uFile2, double time2, char* uFile,     double time);
  void uFileInterpolateCubic     (char* uFile1, double time1, char* uFile2, double time2, char* uFile3,    double time3, char* uFile, double time);
  void uFileInterpolateQuartic   (char* uFile1, double time1, char* uFile2, double time2, char* uFile3,    double time3,
				  char* uFile4, double time4, char* uFile,  double time);

  // Interpolate mFiles
  void mFileInterp               (vector<string> mFileName,   vector<double> mFileTime,   int mFilenumber, double TIME);
  void mFileInterpolateLinear    (char* uFile1, double time1, char* uFile,  double time);
  void mFileInterpolateQuadratic (char* mFile1, double time1, char* mFile2, double time2, char* mFile,     double time);
  void mFileInterpolateCubic     (char* mFile1, double time1, char* mFile2, double time2, char* mFile3,    double time3, char* mFile, double time);
  void mFileInterpolateQuartic   (char* mFile1, double time1, char* mFile2, double time2, char* mFile3,    double time3,
				  char* mFile4, double time4, char* mFile,  double time);

  // Interpolate lFiles
  void lFileInterp               (vector<string> lFileName,   vector<double> lFileTime,   int lFilenumber, double TIME);
  void lFileInterpolateLinear    (char* uFile1, double time1, char* uFile,  double time);
  void lFileInterpolateQuadratic (char* lFile1, double time1, char* lFile2, double time2, char* lFile,     double time);
  void lFileInterpolateCubic     (char* lFile1, double time1, char* lFile2, double time2, char* lFile3,    double time3, char* lFile, double time);
  void lFileInterpolateQuartic   (char* lFile1, double time1, char* lFile2, double time2, char* lFile3,    double time3,
				  char* lFile4, double time4, char* lFile,  double time);

  // Interpolate fFiles
  void fFileInterp               (vector<string> fFileName,   vector<double> fFileTime,   int fFilenumber, double TIME);
  void fFileInterpolateLinear    (char* fFile1, double time1, char* fFile,  double time);
  void fFileInterpolateQuadratic (char* fFile1, double time1, char* fFile2, double time2, char* fFile,     double time);
  void fFileInterpolateCubic     (char* fFile1, double time1, char* fFile2, double time2, char* fFile3,    double time3, char* fFile, double time);
  void fFileInterpolateQuartic   (char* fFile1, double time1, char* fFile2, double time2, char* fFile3,    double time3,
				  char* fFile4, double time4, char* fFile,  double time);
  
  // Interpolate nFiles
  void nFileInterp               (vector<string> nFileName,   vector<double> nFileTime,   int nFilenumber, double TIME);
  void nFileInterpolateLinear    (char* nFile1, double time1, char* nFile,  double time);
  void nFileInterpolateQuadratic (char* nFile1, double time1, char* nFile2, double time2, char* nFile,     double time);
  void nFileInterpolateCubic     (char* nFile1, double time1, char* nFile2, double time2, char* nFile3,    double time3, char* nFile, double time);
  void nFileInterpolateQuartic   (char* nFile1, double time1, char* nFile2, double time2, char* nFile3,    double time3,
				  char* nFile4, double time4, char* nFile,  double time);

  // Find width in PsiN of magnetic island chain
  double GetIslandWidth (int j);
  // Find width in PsiN of vacuum magnetic island chain
  double GetVacuumIslandWidth (int j);
  // Find width in PsiN of vacuum magnetic island chain
  double GetVacuumIslandWidth (int j, double chi);
  // Find limits of magnetic island chains in PsiN
  void GetIslandLimits (int j, double Psi, double& Xminus, double& Xplus);
  // Solve island width equation
  double GetIslandRoot (double c);

  // Open file for reading
  FILE* OpenFiler (char* filename);
  // Open file for writing
  FILE* OpenFilew (char* filename);
  // Open file for appending
  FILE* OpenFilea (char* filename);
  // Call operating system 
  void CallSystem (char* command);
};

#endif //PHASE
