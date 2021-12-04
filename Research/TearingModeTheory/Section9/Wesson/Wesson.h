// Calculate satuarted island widths for Wesson-profile plasma

#ifndef WESSON
#define WESSON

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <blitz/array.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_gamma.h>

using namespace blitz;

// ############
// Class header
// ############
class Wesson
{
 private:
  // ------------------
  // Physics parameters
  // ------------------
  double mpol; // Poloidal mode number
  double ntor; // Toroidal mode number
  double rw;   // Radius of wall (relative to minor radius of plasma)
  double q0;   // Safety-factor at magnetic axis
  double qa;   // Safety factor at edge of plasma 
  
  // ----------------------
  // Integration parameters
  // ----------------------
  double eps;   // Closest distance to magnetic axis
  double delta; // Closest distance to rational surface
  double kmax;  // Integration stopped at k = kmax
  double h0;    // Initial step-length
  double acc;   // Integration accuracy
  double hmin;  // Minimum step-length

  // ----
  // Misc
  // ----
  int count;
  
  // ----------------------
  // Public class functions
  // ----------------------
 public:
  
  Wesson ();              // Constructor
  virtual ~Wesson () {};  // Destructor

  void Solve  ();            // Solve problem

  // -----------------------
  // Private class functions
  // -----------------------
 private:

  // Calculate tearing stability index
  double GetDelta ();
  // Find resonant surface radius
  double Findrs (double qs);
  // Get value fo parameter nu
  double Getnu ();
  // Get current profile
  double GetJ (double r);
  // Get first derivative of current profile
  double GetJp (double r);
  // Get second derivaitive of current profile
  double GetJpp (double r);
  // Get safety-factor profile
  double Getq (double r);
  // Get magnetic shear profile
  double Gets (double r);
  // Get value of parameter alpha
  double Getalpha (double r);
  // Get value of parameter beta
  double Getbeta (double r);

  // Evaluate right-hand sides of differential equations
  void Rhs (double x, Array<double,1>& y, Array<double,1>& dydx);
  // Adaptive step integration routine
  void RK4Adaptive (double& x, Array<double,1>& y, double& h, double& t_err, 
		    double acc, double S, int& rept, int maxrept, 
		    double h_min, double h_max, int flag, int diag, FILE* file);
  // Fixed step integration routine
  void RK4Fixed (double& x, Array<double,1>& y, double h);

  FILE* OpenFile (char* filename);
};

#endif //WESSON
