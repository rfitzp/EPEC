// Neoclassical.h

// ################################################################
// Class to calculate neoclassical flow-damping rates, neoclassical
// phase velocities, and neoclassical resistivities, at rational 
// surfaces in tokamak plasma.

// .....................
// Command line options:
// .....................

// -e INTF     - override INTF value from namelist file
// -p INTP     - override INTP value from namelist file
// -c INTC     - override INTC value from namelist file
// -n NEUTRAL  - override NEUTRAL value from namelist file
// -I IMPURITY - override IMPURITY value from namelist file
// -f FREQ     - override FREQ value from namelist file
// -T NTYPE    - override NTYPE value from namelist file
// -y YN       - override YN value from namelist file
// -N NN       - override NN value from namelist file
// -l LN       - override LN value from namelist file
// -t TIME     - sets experimental time (ms)

// ........
// Outputs:
// ........

// Intermediate data in folder Outputs/Stage3
// Final data passed to program PHASE in file Outputs/nFile

// .........
// Versions:
// .........

// 1.0  - Initial version
// 1.1  - Improved indexing of pFiles and fFiles
// 1.2  - Major rearrangement of input and output files
// 1.3  - Added PsiNk to nFile
// 1.4  - Added cFile and linear interpolation
// 1.5  - Added chie and chin
// 1.6  - Divided normalized layer width by 0.8227
// 1.7  - Output wnl
// 1.8  - Removed 0.8227 from layer width (which actually is not normalized), redefined S_k
// 1.9  - Minor changes for KSTAR data
// 1.10 - Changed operation of FREQ switch
// 1.11 - Added missing extra term associated with chanrge exchange
// 1.12 - Updated fFile input for additonal terms
// 1.13 - Corrected error in fFile interpolation. Improved system calls.
// 1.14 - Renamed Namelist. Adjusted for new fFile format

// ################################################################

#ifndef NEOCLASSICAL
#define NEOCLASSICAL

#define VERSION_MAJOR 1
#define VERSION_MINOR 14

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <blitz/array.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_poly.h>

#include "Field.h"

#define MAXFILENAMELENGTH 500
#define MAXPFILELINELENGTH 500

using namespace blitz;

// Namelist funtion
extern "C" void NameListRead (int* IMPURITY, int* NEUTRAL, int* FREQ, int* INTP, int* INTF, int* INTC, 
			      int* NTYPE, double* NN, double* LN, double* SVN, double* YN, double* EN, double* TIME, double* COULOMB);

// ############
// Class header
// ############
class Neoclassical
{
 private:
  // ++++++++++++++++
  // Class parameters
  // ++++++++++++++++
  
  // ------------------
  // Control parameters
  // ------------------

  // Read from Inputs/Neoclassical.nml
  int    IMPURITY; // Impurity switch. If != 0 then single impurity species included in calculation
  int    NEUTRAL;  // Neutral switch.  If != 0 then majority ion neutrals included in calculation
  int    FREQ;     // Natural frequency switch:
                   //                 If == +1 then use linear frequency derived from given ExB profile
                   //                 If == +2 then use ExB frequency derived from given ExB profile
                   //                 If == +3 then use nonlinear frequency derived from given ExB profile
                   //                 If == -1 then use linear frequency derived from inferred ExB profile
                   //                 If == -2 then use ExB frequency derived from inferred ExB profile
                   //                 If == -3 then use nonlinear derived frequency from inferred ExB profile
  int    INTP;     // If != 0 then use interpolated pFile
  int    INTF;     // If != 0 then use interpolated fFile
  int    INTC;     // If != 0 then use interpolated cFile
  int    NTYPE;    // If == 0 then neutral density distribution exponential. If == 1 then neutral density distribution Lorentzian.
  double NN;       // Flux-surface averaged majority neutral density at plasma boundary (PSI=1) (m^-3)
  double LN;       // Flux-surface averaged majority neutral density decay lengthscale (m)
  double SVN;      // Majority ion/neutral charge exchange rate constant (m^3 /s)
  double YN;       // Majority neutral peaking factor on flux-surfaces
  double EN;       // Ratio of majority neutral to ion temperatures
  double TIME;     // Experimental time (ms)
 
  // ------------------
  // Physical constants
  // ------------------
  double e;         // Magnitude of electron charge
  double epsilon_0; // Electric permittivity of free space
  double mu_0;      // Magnetic permeability of free space
  double m_p;       // Mass of proton
  double m_e;       // Mass of electron
  double COULOMB;   // Coulomb logarithm

  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double xmin;     // Velocity space integrals integrated from x = xmin
  double xmax;     // Velocity space integrals integrated to x = xmax
  double h0;       // Initial step-length
  double acc;      // Integration accuracy
  double hmin;     // Minimum step-length
  double hmax;     // Maximum step-length
  int    maxrept;  // Maximum number of step recalculations

  // .......................
  // Plasma equilibrium data
  // .......................
  
  // From Inputs/fFile
  double          B_0;    // Toroidal magnetic field (T)
  double          R_0;    // Major radius (m)
  double          a;      // Minor radius (m)
  double          h_a;    // a/R_0
  int             NPSI;   // Number of equilibrium data points
  Array<double,1> psi;    // Normalized poloidal magnetic flux
  Array<double,1> rr;     // Normalized flux-surface minor radius
  Array<double,1> dpsidr; // Normalized dpsi/dr

  // -------------------
  // Plasma profile data
  // -------------------

  // Read from Inputs/pfile
  Field ne;   // Electron number density (m^-3)
  Field Te;   // Electron temperature (J)
  Field ni;   // Majority ion number density (m^-3)
  Field Ti;   // Majority ion temperature (J)
  Field nI;   // Impurity ion number density (m^-3)
  Field nb;   // Fast majority ion number density (m^-3)
  Field wE;   // ExB frequency (rad/s)
  Field wt;   // Impurity ion toroidal angular frequency (rad/s)
  Field NZA;  // Ion data

  int    NI;  // Charge number of majority ions
  double ZI;  // Ionization number of majority ions
  double AI;  // Mass number of majority ions
  int    NII; // Charge number of impurity ions
  double ZII; // Ionization number of impurity ions
  double AII; // Mass number of impurity ions

  // Read from Inputs/cFile
  Field Chip; // Perpendicular momentum diffusivity
  Field Chie; // Perpendicular energy   diffusivity
  Field Chin; // Perpendicular particle diffusivity

  // Profile data interpolated onto equilibrium grid
  Array<double,1> n_e;    // Electron number density (m^-3)
  Array<double,1> dn_edr; // Electron number density gradient (m^-4)
  Array<double,1> T_e;    // Electron temperature  (J)
  Array<double,1> dT_edr; // Electron temperature gradient (J m^-1)
  Array<double,1> n_i;    // Majority ion number density (m^-3)
  Array<double,1> dn_idr; // Majority ion number density gradient  (m^-4)
  Array<double,1> T_i;    // Majority ion temperature (J)
  Array<double,1> dT_idr; // Majority ion temperature gradient (J m^-1)
  Array<double,1> n_b;    // Fast majority ion number density (m^-3)
  Array<double,1> n_I;    // Impurity ion number density (m^-3)
  Array<double,1> dn_Idr; // Impurity ion number density gradient (m^-4)
  Array<double,1> T_I;    // Impurity ion temperature (J)
  Array<double,1> dT_Idr; // Impurity ion temperature gradient (J m^-1)
  Array<double,1> n_n;    // Neutral number density (m^-3)
  Array<double,1> w_E;    // ExB frequency (rad/s)
  Array<double,1> w_t;    // Impurity ion toroidal angular frequency (rad/s)
  Array<double,1> Quasi;  // Quasi-nuetrality check
  Array<double,1> Z_eff;  // Effective ion charge number
  Array<double,1> alpha;  // Impurity strength parameter
  Array<double,1> chip;   // Perpendicular momentum diffusivity (m^2 s^-1)
  Array<double,1> chie;   // Perpendicular energy   diffusivity (m^2 s^-1)
  Array<double,1> chin;   // Perpendicular particle diffusivity (m^2 s^-1)

  // ---------------------------
  // Rational surface parameters
  // ---------------------------

  // Read from Inputs/fFile
  int             ntor;   // Toroidal mode number
  int             nres;   // Number of resonant surfaces
  Array<int,1>    mk;     // Poloidal mode numbers
  Array<double,1> rk;     // Normalized minor radii
  Array<double,1> PsiNk;  // Normalized polidal fluxes at rational surface
  Array<double,1> qk;     // Safety factors
  Array<double,1> sk;     // Magnetic shears
  Array<double,1> gk;     // g values
  Array<double,1> gmk;    // gamma values
  Array<double,1> Ktk;    // Kt values
  Array<double,1> Kastk;  // Kast values
  Array<double,1> Kthek;  // Ktheta values
  Array<double,1> fck;    // Fractions of circulating particles
  Array<double,1> akk;    // Metric elements
  Array<double,1> dPsidr; // dPsiN/dr

  // Derived from profiles
  double rho0;                // Central mass density
  double tau_A;               // Central Alfven time
  
  Array<double,1> nek;        // Electron number densities
  Array<double,1> dnedrk;     // Electron number density gradients
  Array<double,1> Tek;        // Electron temperatures
  Array<double,1> dTedrk;     // Electron temperature gradients
  Array<double,1> nik;        // Majority ion number densities
  Array<double,1> dnidrk;     // Majority ion number density gradients
  Array<double,1> Tik;        // Majority ion temperatures
  Array<double,1> dTidrk;     // Majority ion temperature gradients
  Array<double,1> nbk;        // Fast majority ion numbers densities
  Array<double,1> nIk;        // Impurity ion number densities
  Array<double,1> dnIdrk;     // Impurity ion number density gradients
  Array<double,1> TIk;        // Impurity ion temperatures
  Array<double,1> dTIdrk;     // Impurity ion temperature gradients
  Array<double,1> wEk;        // ExB frequencies
  Array<double,1> wtk;        // Impurity ion toroidal rotation frequencies
  Array<double,1> Zeffk;      // Effective ion charge numbers
  Array<double,1> alphak;     // Impurity strength parameters
  Array<double,1> rhok;       // Relative mass densities
  Array<double,1> NNk;        // Flux-surfaced-averaged majority neutral number densities
  Array<double,1> chipk;      // Perpendicular momentum diffusivities
  Array<double,1> chiek;      // Perpendicular energy   diffusivities
  Array<double,1> chink;      // Perpendicular particle diffusivities

  Array<double,1> v_T_ek;     // Electron thermal velocities
  Array<double,1> v_T_ik;     // Majority ion thermal velocities
  Array<double,1> v_T_Ik;     // Impurity ion thermal velocities

  Array<double,1> omega_t_ek; // Electron transit frequencies
  Array<double,1> omega_t_ik; // Majority ion transit frequencies
  Array<double,1> omega_t_Ik; // Impurity ion transit frequencies

  Array<double,1> nu_eek;     // Electron collision frequencies  
  Array<double,1> nu_iik;     // Majority ion collision frequencies
  Array<double,1> nu_IIk;     // Impurity ion collision frequencies

  Array<double,1> WcritTk;    // Critical island width for temperature flattening
  Array<double,1> Wcritnk;    // Critical island width for density flattening

  Array<double,1> eta_ek;     // Relative electron temperature gradients
  Array<double,1> eta_ik;     // Relative majority ion temperature gradients
  Array<double,1> eta_Ik;     // Relative impurity ion temperature gradients

  Array<double,1> w_ast_ek;   // Electron diamagnetic frequencies
  Array<double,1> w_ast_ik;   // Majority ion diamagnetic frequencies
  Array<double,1> w_ast_Ik;   // Impurity ion diamagnetic frequencies
  Array<double,1> w_nc_ik;    // Majority ion neoclassical frequencies
  Array<double,1> w_nc_Ik;    // Impurity ion neoclassical frequencies
  Array<double,1> w_E_Ik;     // ExB frequencies inferred from toroidal impurity ion rotation frequency

  Array<double,1> rho_sk;     // Ion sound radii
  
  Array<double,1> tau_Hk;     // Hydromagnetic timescales
  Array<double,1> tau_Rk;     // Classical resistive timescales
  Array<double,1> tau_Mk;     // Momentum confinement timescales
  Array<double,1> tau_thk;    // Poloidal flow damping timescales

  Array<double,1> gt;         // Fraction of trapped particles
  Array<double,1> nu_P_e;     // Electron bananna/plateau collisionality parameter
  Array<double,1> nu_P_i;     // Majority ion bananna/plateau collisionality  parameter
  Array<double,1> nu_P_I;     // Impurity ion bananna/plateau collisionality parameter
  Array<double,1> nu_PS_e;    // Electron plateau/Pfirsch-Schluter collisionality parameter
  Array<double,1> nu_PS_i;    // Majority ion plateau/Pfirsch-Schluter collisionality parameter
  Array<double,1> nu_PS_I;    // Impurity ion plateau/Pfirsch-Schluter collisionality parameter
  Array<double,1> x_iI;       // Ratios of ion thermal speeds
  Array<double,1> x_Ii;       // Inverse ratios of ion thermal speeds
  
  // ------------------------
  // Neoclassical viscosities
  // ------------------------
  Array<double,1> mu_00_i;  // Majority ion neoclassical viscosities
  Array<double,1> mu_01_i;  // Majority ion Neoclassical viscosities
  Array<double,1> mu_11_i;  // Majority ion Neoclassical viscosities
  Array<double,1> mu_00_I;  // Impurity ion neoclassical viscosities
  Array<double,1> mu_01_I;  // Impurity ion neoclassical viscosities
  Array<double,1> mu_11_I;  // Impurity ion neoclassical viscosities
  Array<double,1> mu_00_e;  // Electron neoclassical viscosities
  Array<double,1> mu_01_e;  // Electron neoclassical viscosities
  Array<double,1> mu_11_e;  // Electron neoclassical viscosities

  // -----------------------
  // Neoclassical parameters
  // -----------------------
  Array<double,1> L_ii_00;  // Neoclassical ion flow parameter
  Array<double,1> L_ii_01;  // Neoclassical ion flow parameter
  Array<double,1> L_iI_00;  // Neoclassical ion flow parameter
  Array<double,1> L_iI_01;  // Neoclassical ion flow parameter
  Array<double,1> L_Ii_00;  // Neoclassical ion flow parameter
  Array<double,1> L_Ii_01;  // Neoclassical ion flow parameter
  Array<double,1> L_II_00;  // Neoclassical ion flow parameter
  Array<double,1> L_II_01;  // Neoclassical ion flow parameter
  Array<double,1> G_ii_00;  // Neoclassical ion flow parameter
  Array<double,1> Q_00;     // Neoclassical electron flow parameter

  // ........................
  // Neoclassical frequencies
  // ........................
  Array<double,1> w_linear;     // Neoclassical frequencies from linear theory
  Array<double,1> w_nonlinear;  // Neoclassical frequencies from nonlinear theory
  Array<double,1> w_EB;         // Neoclassical frequencies assuming convection by ExB fluid
  Array<double,1> w_actual;     // Actual neoclassical frequencies
  Array<double,1> w_fac;        // Degree of island propagation in ion diamagnetic direction
  
  // ----
  // Misc
  // ----
  int jj, count, swit;
  
  // ----------------------
  // Public class functions
  // ----------------------
 public:
  
  Neoclassical ();             // Constructor
  virtual ~Neoclassical () {}; // Destructor

  // Solve problem
  void Solve  (int _NEUTRAL, int _IMPURITY, int _FREQ, int _INTP, int _INTF,
	       int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME);            
 
  // -----------------------
  // Private class functions
  // -----------------------
 private:

  // Read discharge parameters
  void Read_Parameters (int _NEUTRAL, int _IMPURITY, int _FREQ, int _INTP, int _INTF,
			int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME);
  // Read equilibrium data
  void Read_Equilibrium ();
  // Read profile data
  void Read_Profiles ();
  // Read pFile
  void pFileRead ();
  // Read cFile
  void cFileRead ();
  // Calculate derived quantities at rational surfaces
  void Get_Derived ();
  // Calculate neoclassical viscosities at sational surfacea
  void Get_Viscosities ();
  // Calculate neoclassical parameters at sational surfacea
  void Get_Parameters ();
  // Calculate neoclassical freqeuncies at rational surfaces
  void Get_Frequencies ();
  // Calculate normalized quantities at rational surfaces
  void Get_Normalized ();

  // 1D interpolation function with nonuniform grid
  double Interpolate             (int I, Array<double,1> X, Array<double,1> Y, double x, int order);
  double InterpolateCubic        (Array<double,1> X, Array<double,1> Y, double x, int i0, int i1, int i2,         int order);
  double InterpolateQuartic      (Array<double,1> X, Array<double,1> Y, double x, int i0, int i1, int i2, int i3, int order);
  double InterpolateField        (Field& F, double x, int order);
  double InterpolateFieldCubic   (Field& F, double x, int i0, int i1, int i2,         int order);
  double InterpolateFieldQuartic (Field& F, double x, int i0, int i1, int i2, int i3, int order);

  // Interpolate pFiles
  void pFileInterp               (vector<string> pFileName, vector<double> pFileTime, int pFileNumber, double time);
  void pFileInterpolateLinear    (char* pFile1, double time1, char* pFile,  double time);
  void pFileInterpolateQuadratic (char* pFile1, double time1, char* pFile2, double time2, char* pFile, double time);
  void pFileInterpolateCubic     (char* pFile1, double time1, char* pFile2, double time2, char* pFile3, double time3, char* pFile,  double time);
  void pFileInterpolateQuartic   (char* pFile1, double time1, char* pFile2, double time2, char* pFile3, double time3, char* pFile4, double time4,
				  char* pFile,  double time);

  // Interpolate cFiles
  void cFileInterp               (vector<string> cFileName, vector<double> cFileTime, int cFileNumber, double time);
  void cFileInterpolateLinear    (char* cFile1, double time1, char* cFile,  double time);
  void cFileInterpolateQuadratic (char* cFile1, double time1, char* cFile2, double time2, char* cFile, double time);
  void cFileInterpolateCubic     (char* cFile1, double time1, char* cFile2, double time2, char* cFile3, double time3, char* cFile,  double time);
  void cFileInterpolateQuartic   (char* cFile1, double time1, char* cFile2, double time2, char* cFile3, double time3, char* cFile4, double time4,
				  char* cFile,  double time);
  
  // Interpolate Fields
  void FieldInterpolateLinear    (Field& Field1, Field& Field,  double weight1);
  void FieldInterpolateQuadratic (Field& Field1, Field& Field2, Field& Field,
				  double weight1, double weight2);
  void FieldInterpolateCubic     (Field& Field1, Field& Field2, Field& Field3, Field& Field,
				  double weight1, double weight2, double weight3);
  void FieldInterpolateQuartic   (Field& Field1, Field& Field2, Field& Field3, Field& Field4, Field& Field,
				  double weight1, double weight2, double weight3, double weight4);

  // Interpolate fFiles
  void fFileInterp               (vector<string> fFileName, vector<double> fFileTime, int fFileNumber, double time);
  void fFileInterpolateLinear    (char* fFile1, double time1, char* fFile,  double time);
  void fFileInterpolateQuadratic (char* fFile1, double time1, char* fFile2, double time2, char* fFile, double time);
  void fFileInterpolateCubic     (char* fFile1, double time1, char* fFile2, double time2, char* fFile3, double time3, char* fFile, double time);
  void fFileInterpolateQuartic   (char* fFile1, double time1, char* fFile2, double time2, char* fFile3, double time3,
				  char* fFile4, double time4, char* fFile,  double time);
 
   // Chandrasekhar function
  double psi_fun (double x);
  // Derivative of Chandrasekhar function
  double psi_fun_p (double x);
   
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, Array<double,1>& y, Array<double,1>& dydx);
  // Adaptive-step integration routine
  void RK4Adaptive (double& x, Array<double,1>& y, double& h, double& t_err, 
		    double acc, double S, int& rept, int maxrept, 
		    double h_min, double h_max, int flag, int diag, FILE* file);
  // Fixed step integration routine
  void RK4Fixed (double& x, Array<double,1>& y, double h);

  // Open file for reading
  FILE* OpenFiler (char* filename);
  // Open file for writing
  FILE* OpenFilew (char* filename);
  // Open file for appending
  FILE* OpenFilea (char* filename);
  // Call operating system 
  void CallSystem (char* command);
};

#endif //NEOCLASSICAL
