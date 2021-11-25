// Calculate integrals in Section 8

#ifndef INTEGRALS
#define INTEGRALS

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
class Integrals
{
 private:
  // ----------------------
  // Integration parameters
  // ----------------------
  double eps;   // Integration started at k = 1 + eps
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
  
  Integrals ();              // Constructor
  virtual ~Integrals () {};  // Destructor

  void Solve  ();            // Solve problem

  // -----------------------
  // Private class functions
  // -----------------------
 private:

  // Useful functions
  double Ak  (double k);
  double Bk  (double k);
  double Ck  (double k);
  double Dk  (double k);
  double Ek  (double k);

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

#endif //INTEGRALS
