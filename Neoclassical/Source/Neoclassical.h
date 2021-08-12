// Neoclassical.h

// ################################################################
// Class to calculate neoclassical flow-damping rates, neoclassical
// phase velocities, and neoclassical resistivities, at rational 
// surfaces in tokamak plasma.

// .....................
// Command line options:
// .....................

// -c INTC     - override INTC value from namelist file
// -e INTF     - override INTF value from namelist file
// -f EXB      - override EXB value from namelist file
// -h          - list options
// -l LN       - override LN value from namelist file
// -o          - flag to select OMFIT mode
// -n NEUTRAL  - override NEUTRAL value from namelist file
// -p INTP     - override INTP value from namelist file
// -t TIME     - sets experimental time (ms)
// -y YN       - override YN value from namelist file
// -C CATS     - override CATS value from namelist file
// -I IMPURITY - override IMPURITY value from namelist file
// -N NN       - override NN value from namelist file
// -T NTYPE    - override NTYPE value from namelist file

// ...................
// Inputs and outputs:
// ...................
// Calculation control parameters in namelist file Inputs/Neoclasssical.nml

// FLUX data in Inputs/fFile or Inputs/fFiles
// Profile data in Inputs/pFile or Inputs pFiles
// Transport data in Inputs/cFile or Inputs/cFiles

// Intermediate data in folder Outputs/Stage3
// Final data passed to program PHASE in file Outputs/nFile or Outputs/nFiles

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
// 1.15 - Added extra information to nFile
// 1.16 - Added P0 to nFile
// 1.17 - Added higher order derivatives of profiles
// 1.18 - Added L_Ii_00 correction
// 1.19 - Added Chii. Changed FREQ flag to EXB flag.
// 1.20 - Added CATS flag
// 1.21 - Use more accurate expression for neoclassical flow damping time
// 1.22 - Added more accurate calculation of linear layer width
// 1.23 - Added C1 and C2
// 1.24 - Added change exchange to angular velocity evolution equations
// 1.25 - Added calculation of DB and DR
// 1.26 - Added calculation of alpha_b(e,i), alpha_c, and alpha_p
// 1.27 - Added NETCDF output. Adapted ti run with OMFIT.

// ################################################################

#ifndef NEOCLASSICAL
#define NEOCLASSICAL

#define VERSION_MAJOR 1
#define VERSION_MINOR 27
#define MAXFILENAMELENGTH 500
#define MAXPFILELINELENGTH 500

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <blitz/array.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_poly.h>
#include "Field.h"
#include <netcdf.h>

using namespace blitz;

// Namelist function
extern "C" void NameListRead (int* IMPURITY, int* NEUTRAL, int* EXB, int* INTP, int* INTF, int* INTC, 
			      int* NTYPE, double* NN, double* LN, double* SVN, double* YN, double* EN,
			      double* TIME, double* COULOMB, int* NSMOOTH, int *CATS, double* TAUMIN);

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

  // Read from command line
  int OMFIT;       // Flag to select OMFIT mode

  // Read from Inputs/Neoclassical.nml
  double COULOMB;  // Coulomb logarithm
  int    NSMOOTH;  // Number of smoothing cycles for higher devivatives of profiles

  int    INTP;     // If != 0 then use interpolated pFile
  int    INTF;     // If != 0 then use interpolated fFile
  int    INTC;     // If != 0 then use interpolated cFile
  int    CATS;     // If != 0 then use only linear interpolation for cFiles

  int    EXB;      //  If == 0 then use ExB frequency profile from pFile
                   //  If == 1 then use ExB frequency profile derived from pFile toroidal velocity profile and neoclassical theory  

  int    IMPURITY; // Impurity switch. If != 0 then single impurity species included in calculation

  int    NEUTRAL;  // Neutral switch.  If != 0 then majority ion neutrals included in calculation
  int    NTYPE;    // If == 0 then neutral density distribution exponential. If == 1 then neutral density distribution Lorentzian.
  double NN;       // Flux-surface averaged majority neutral density at plasma boundary (PSI=1) (m^-3)
  double LN;       // Flux-surface averaged majority neutral density decay lengthscale (m)
  double SVN;      // Majority ion/neutral charge exchange rate constant (m^3 /s)
  double YN;       // Majority neutral peaking factor on flux-surfaces
  double EN;       // Ratio of majority neutral to ion temperatures

  double TIME;     // Experimental time (ms)

  double TAUMIN;   // Minimum allowed value of tau (i.e., ratio of electon to ion diamagnetic frequencies)
                   //  1+tau must be positive otherwise layer width calculation fails.

  // ------------------
  // Physical constants
  // ------------------
  double e;         // Magnitude of electron charge (SI)
  double epsilon_0; // Electric permittivity of free space (SI)
  double mu_0;      // Magnetic permeability of free space (SI)
  double m_p;       // Mass of proton (SI)
  double m_e;       // Mass of electron (SI)

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
  int    flag;     // Flag for selecting right-hand sides of differential equations
  int    jres;     // Flag for selecting resonant surface

  // .......................
  // Plasma equilibrium data
  // .......................
  
  // From Inputs/fFile
  double          B_0;    // Toroidal magnetic field (T)
  double          R_0;    // Major radius (m)
  double          a;      // Minor radius (m)
  double          h_a;    // a /R_0
  int             NPSI;   // Number of equilibrium data points
  Array<double,1> psi;    // Normalized poloidal magnetic flux (0 at axis, 1 at lcfs)
  Array<double,1> rr;     // Flux-surface minor radius / a
  Array<double,1> dpsidr; // a dpsi/dr
  double          PSIRAT; // Value of PsiN beyond which resonant surfaces ignored

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
  Field Chip; // Perpendicular momentum diffusivity (m^2/s)
  Field Chie; // Perpendicular electron energy diffusivity (m^2/s)
  Field Chin; // Perpendicular particle diffusivity (m^2/s)
  Field Chii; // Perpendicular ion energy diffusivity (m^2/s)

  // Profile data interpolated onto equilibrium grid
  Array<double,1> n_e;    // Electron number density (m^-3)
  Array<double,1> dn_edr; // Electron number density gradient (in r) (m^-4)
  Array<double,1> T_e;    // Electron temperature (J)
  Array<double,1> dT_edr; // Electron temperature gradient (in r) (J m^-1)
  Array<double,1> n_i;    // Majority ion number density (m^-3)
  Array<double,1> dn_idr; // Majority ion number density gradient (in r) (m^-4)
  Array<double,1> T_i;    // Majority ion temperature (J)
  Array<double,1> dT_idr; // Majority ion temperature gradient (in r) (J m^-1)
  Array<double,1> n_b;    // Fast majority ion number density (m^-3)
  Array<double,1> n_I;    // Impurity ion number density (m^-3)
  Array<double,1> dn_Idr; // Impurity ion number density gradient (m^-4)
  Array<double,1> T_I;    // Impurity ion temperature (J)
  Array<double,1> dT_Idr; // Impurity ion temperature gradient (in r) (J m^-1)
  Array<double,1> n_n;    // Neutral number density (m^-3)
  Array<double,1> w_E;    // ExB frequency (rad/s)
  Array<double,1> w_t;    // Impurity ion toroidal angular frequency (rad/s)
  Array<double,1> Quasi;  // Quasi-nuetrality check
  Array<double,1> Z_eff;  // Effective ion charge number
  Array<double,1> alpha;  // Impurity strength parameter
  Array<double,1> chip;   // Perpendicular momentum diffusivity (m^2 s^-1)
  Array<double,1> chie;   // Perpendicular electron energy diffusivity (m^2 s^-1)
  Array<double,1> chin;   // Perpendicular particle diffusivity (m^2 s^-1)
  Array<double,1> chii;   // Perpendicular ion energy diffusivity (m^2 s^-1)

  Array<double,1> dn_edP1; // Electron number density 1st derivative in PsiN (m^-3)
  Array<double,1> dT_edP1; // Electron temperature 1st derivative in PsiN (J)
  Array<double,1> dn_idP1; // Majority ion number density 1st derivative in PsiN (m^-3)
  Array<double,1> dT_idP1; // Majority ion temperature 1st derivative in PsiN (J)
  Array<double,1> dn_edP2; // Electron number density 2nd derivative in PsiN (m^-3)
  Array<double,1> dT_edP2; // Electron temperature 2nd derivative in PsiN (J)
  Array<double,1> dn_idP2; // Majority ion number density 2nd derivative in PsiN (m^-3)
  Array<double,1> dT_idP2; // Majority ion temperature 2nd derivative in PsiN (J)
  Array<double,1> dn_edP3; // Electron number density 3rd derivative in PsiN (m^-3)
  Array<double,1> dT_edP3; // Electron temperature 3rd derivative in PsiN (J)
  Array<double,1> dn_idP3; // Majority ion number density 3rd derivative in PsiN (m^-3)
  Array<double,1> dT_idP3; // Majority ion temperature 3rd derivative in PsiN (J)

  Array<double,1> Factor1;  // First transport factor
  Array<double,1> Factor2;  // Second transport factor
  Array<double,1> Factor3;  // Third transport factor
  Array<double,1> Factor4;  // Fourth transport factor
  Array<double,1> Factor5;  // Fifth transport factor
  Array<double,1> Factor6;  // Sixth transport factor
  Array<double,1> Factor7;  // Seventh transport factor
  Array<double,1> Factor8;  // Eight transport factor
  Array<double,1> Factor9;  // Ninth transport factor
  Array<double,1> Factor10; // Tenth transport factor
  Array<double,1> Factor11; // Eleventh transport factor
  Array<double,1> Factor12; // Twelth transport factor
  
  // ---------------------------
  // Rational surface parameters
  // ---------------------------

  // Read from Inputs/fFile
  int             ntor;   // Toroidal mode number
  int             nres;   // Number of resonant surfaces
  Array<int,1>    mk;     // Poloidal mode numbers
  Array<double,1> rk;     // Minor radii /a
  Array<double,1> PsiNk;  // Normalized poloidal fluxes 
  Array<double,1> qk;     // Safety factors
  Array<double,1> sk;     // Magnetic shears (dlnq/dlnr)
  Array<double,1> gk;     // g values
  Array<double,1> gmk;    // gamma values
  Array<double,1> Ktk;    // Kt values
  Array<double,1> Kastk;  // Kast values
  Array<double,1> Kthek;  // Ktheta values
  Array<double,1> fck;    // Fractions of circulating particles
  Array<double,1> akk;    // Metric elements
  Array<double,1> dPsidr; // R_0 dPsiN/dr
  Array<double,1> A2;     // A2
  Array<double,1> q_hat;  // q_hat
  Array<double,1> C1;     // C1
  Array<double,1> C2;     // C2
  Array<double,1> DR;     // GGJ stability parameter

  // Derived from profiles
  double rho0;                // Central mass density (kg/m^-3)
  double tau_A;               // Central Alfven time (s)
  double P0;                  // Central (thermal) pressure (Pa)
  
  Array<double,1> nek;        // Electron number densities (m^-3)
  Array<double,1> dnedrk;     // Electron number density gradients (in r) (m^-4)
  Array<double,1> Tek;        // Electron temperatures (J)
  Array<double,1> dTedrk;     // Electron temperature gradients (in r) (J m^-1)
  Array<double,1> nik;        // Majority ion number densities (m^-3)
  Array<double,1> dnidrk;     // Majority ion number density gradients (in r) (m^-4)
  Array<double,1> Tik;        // Majority ion temperatures (J)
  Array<double,1> dTidrk;     // Majority ion temperature gradients (in r) (J m^-1)
  Array<double,1> nbk;        // Fast majority ion numbers densities (m^-3)
  Array<double,1> nIk;        // Impurity ion number densities (m^-3)
  Array<double,1> dnIdrk;     // Impurity ion number density gradients (in r) (m^-4)
  Array<double,1> TIk;        // Impurity ion temperatures (J)
  Array<double,1> dTIdrk;     // Impurity ion temperature gradients (in r) (J m^-1)
  Array<double,1> wEk;        // ExB frequencies (rad/s)
  Array<double,1> wtk;        // Impurity ion toroidal rotation frequencies (rad/s)
  Array<double,1> Zeffk;      // Effective ion charge numbers
  Array<double,1> Zeffik;     // Effective ion charge numbers
  Array<double,1> ZeffIk;     // Effective ion charge number
  Array<double,1> alphak;     // Impurity strength parameters
  Array<double,1> rhok;       // Relative mass densities
  Array<double,1> NNk;        // Flux-surfaced-averaged majority neutral number densities (m^-3)
  Array<double,1> chipk;      // Perpendicular momentum diffusivities (m^2/s)
  Array<double,1> chiek;      // Perpendicular electron energy diffusivities (m^2/s)
  Array<double,1> chink;      // Perpendicular particle diffusivities (m^2/s)
  Array<double,1> chiik;      // Perpendicular ion energy diffusivities (m^2/s)
  Array<double,1> alpbek;     // Bootstrap current parameter
  Array<double,1> alpbik;     // Bootstrap current parameter
  Array<double,1> alpck;      // Curvature parameter
  Array<double,1> alppk;      // Polarization current parameter
  Array<double,1> rhothek;    // Electron poloidal gyroradius (m)
  Array<double,1> rhothik;    // Majority ion poloidal gyroradius (m)

  Array<double,1> dnedP1k;    // Electron number density 1st derivative wrt PsiN (m^-3)
  Array<double,1> dTedP1k;    // Electron temperature 1st derivative wrt PsiN (J)
  Array<double,1> dnidP1k;    // Ion number density 1st derivative wrt PsiN (m^-3)
  Array<double,1> dTidP1k;    // Ion temperature 1st derivative wrt PsiN (J)
  Array<double,1> dnedP2k;    // Electron number density 2nd derivative wrt PsiN (m^-3)
  Array<double,1> dTedP2k;    // Electron temperature 2nd derivative wrt PsiN (J)
  Array<double,1> dnidP2k;    // Ion number density 2nd derivative wrt PsiN (m^-3)
  Array<double,1> dTidP2k;    // Ion temperature 2nd derivative wrt PsiN (J)
  Array<double,1> dnedP3k;    // Electron number density 3rd derivative wrt PsiN (m^-3)
  Array<double,1> dTedP3k;    // Electron temperature 3rd derivative wrt PsiN (J)
  Array<double,1> dnidP3k;    // Ion number density 3rd derivative wrt PsiN (m^-3)
  Array<double,1> dTidP3k;    // Ion temperature 3rd derivative wrt PsiN (J)

  Array<double,1> v_T_ek;     // Electron thermal velocities (m/s)
  Array<double,1> v_T_ik;     // Majority ion thermal velocities (m/s)
  Array<double,1> v_T_Ik;     // Impurity ion thermal velocities (m/s)

  Array<double,1> omega_t_ek; // Electron transit frequencies (s^-1)
  Array<double,1> omega_t_ik; // Majority ion transit frequencies (s^-1)
  Array<double,1> omega_t_Ik; // Impurity ion transit frequencies (s^-1)

  Array<double,1> nu_eek;     // Electron collision frequencies (s^-1)
  Array<double,1> nu_iik;     // Majority ion collision frequencies (s^-1)
  Array<double,1> nu_IIk;     // Impurity ion collision frequencies (s^-1)

  Array<double,1> WcritTek;   // Critical island width for electron temperature flattening (in r) (m)
  Array<double,1> WcritTik;   // Critical island width for ion temperature flattening (in r) (m)
  Array<double,1> Wcritnek;   // Critical island width for density flattening (in r) (m)

  Array<double,1> eta_ek;     // Relative electron temperature gradients
  Array<double,1> eta_ik;     // Relative majority ion temperature gradients
  Array<double,1> eta_Ik;     // Relative impurity ion temperature gradients

  Array<double,1> w_ast_ek;   // Electron diamagnetic frequencies (rad/s)
  Array<double,1> w_ast_ik;   // Majority ion diamagnetic frequencies (rad/s)
  Array<double,1> w_ast_Ik;   // Impurity ion diamagnetic frequencies (rad/s)
  Array<double,1> w_nc_ik;    // Majority ion neoclassical frequencies (rad/s)
  Array<double,1> w_nc_Ik;    // Impurity ion neoclassical frequencies (rad/s)
  Array<double,1> w_nc_eek;   // Electron-electron neoclassical frequencies (rad/s)
  Array<double,1> w_nc_eik;   // Electron-ion neocalssical frequencies (rad/s)
  Array<double,1> w_E_Ik;     // ExB frequencies inferred from toroidal impurity ion rotation frequency (rad/s)
  Array<double,1> w_betak;    // Bootstrap frequencies (rad/s)
  Array<double,1> w_Omegk;    // Polarization frequencies (rad/s)

  Array<double,1> rho_sk;     // Ion sound radii (m)
  
  Array<double,1> tau_Hk;     // Hydromagnetic timescales (s)
  Array<double,1> tau_Rk;     // Classical resistive timescales (s)
  Array<double,1> tau_Mk;     // Momentum confinement timescales (s)
  Array<double,1> tau_thk;    // Poloidal flow damping timescales (s)
  Array<double,1> tau_cxk;    // Charge exchange damping timescales (s)

  Array<double,1> Sk;         // Lundquist numbers   
  Array<double,1> tauk;       // Ratio of diamagnetic frequencies
  Array<double,1> PEk;        // Perpendicular particle/energy transport parameter
  Array<double,1> PMk;        // Perpendicular momentum transport parameter
  Array<double,1> Dk;         // Semi-collisional parameter
  Array<double,1> QEk;        // ExB frequency parameter
  Array<double,1> Qek;        // Electron diamagnetic frequency parameter
  Array<double,1> Qik;        // Ion diamagnetic frequency parameter
  Array<double,1> delk;       // Linear layer width (m)
  
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
  Array<double,1> G_ii_00;  // Neoclassical majority ion flow parameter
  Array<double,1> L_ii_00;  // Neoclassical majority ion flow parameter
  Array<double,1> L_ii_01;  // Neoclassical majority ion flow parameter
  Array<double,1> L_iI_00;  // Neoclassical majority ion flow parameter
  Array<double,1> L_iI_01;  // Neoclassical majority ion flow parameter

  Array<double,1> G_Ii_00;  // Neoclassical impurity ion flow parameter
  Array<double,1> L_Ii_00;  // Neoclassical impurity ion flow parameter
  Array<double,1> L_Ii_01;  // Neoclassical impurity ion flow parameter
  Array<double,1> L_II_00;  // Neoclassical impurity ion flow parameter
  Array<double,1> L_II_01;  // Neoclassical impurity ion flow parameter

  Array<double,1> G_ei_00;  // Neoclassical electron flow parameter
  Array<double,1> L_ee_00;  // Neoclassical electron flow parameter
  Array<double,1> L_ee_01;  // Neoclassical electron flow parameter
  Array<double,1> L_ei_00;  // Neoclassical electron flow parameter
  Array<double,1> L_ei_01;  // Neoclassical electron flow parameter
  Array<double,1> L_eI_00;  // Neoclassical electron flow parameter
  Array<double,1> L_eI_01;  // Neoclassical electron flow parameter
  Array<double,1> Q_00;     // Neoclassical electron flow parameter

  // ...................
  // Natural frequencies
  // ...................
  Array<double,1> w_linear;     // Natural frequencies from linear theory (rad/s)
  Array<double,1> w_nonlinear;  // Natural frequencies from nonlinear theory (rad/s)
  Array<double,1> w_EB;         // Natural frequencies assuming convection by ExB fluid (rad/s)
  
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
  void Solve  (int _NEUTRAL, int _IMPURITY, int _EXB, int _INTP, int _INTF,
	       int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME, int _CATS, int _OMFIT);            
 
  // -----------------------
  // Private class functions
  // -----------------------
 private:

  // Read discharge parameters
  void Read_Parameters (int _NEUTRAL, int _IMPURITY, int _EXB, int _INTP, int _INTF,
			int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME, int _CATS, int _OMFIT);
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
  // Calculate linear layer widths at rational surfaces
  void Get_LayerWidths ();
  // Calculate normalized quantities at rational surfaces
  void Get_Normalized ();

  // Write Stage2 NETCDF file
  void WriteStage2Netcdfc ();

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

  // Savitsky-Gorlay smoothing routine
  void Smoothing (int N, Array<double,1> y);

  // Matrix multiplication routine
  void Matrix_Mult (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i);
  // Matrix addition routine
  void Matrix_Add (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i);
  // Matrix subtraction routine
  void Matrix_Sub (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i);

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
