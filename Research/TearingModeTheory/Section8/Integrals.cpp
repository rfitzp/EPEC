// Integrals.cpp

#include "Integrals.h"

// ###########
// Constructor
// ###########
Integrals::Integrals ()
{
  // --------------------------------------
  // Set default values of class parameters
  // --------------------------------------

  kmax  = 400.;
  eps   = 1.e-8;
  h0    = 1.e-6;
  acc   = 1.e-13;
  hmin  = 1.e-10;
}

// #########################
// Function to solve problem
// #########################
void Integrals::Solve ()
{
  // ....................
  // Calculate  integrals
  // ....................

  double k, h, t_err;
  int rept; count = 0;
  Array<double,1> y(4);

  k    = eps;
  y(0) = 0.;
  y(1) = 0.;
  y(2) = 0.;
  y(3) = 0.;
  h    = h0;

  do
    {
      RK4Adaptive (k, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-1, 2, 0, NULL);
    }
  while (k < 1. - eps);

  k = 1. + eps;
  h = h0; 

  do
    {
      RK4Adaptive (k, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-1, 2, 0, NULL);
    }
  while (k < kmax);

  double I1 = 4.             * y(1);
  double I2 = 4. * (1./y(0))      * y(2);
  double I3 = 16. * (1./y(0)/y(0)) * y(3);

  printf ("\nI1 = %11.4e  I2 = %11.4e  I3 = %11.4e\n\n", I1, I2, I3);
}

// ##############################################################
// Function to evaluate right-hand side of differential equations
// ##############################################################
void Integrals::Rhs (double k, Array<double,1>& y, Array<double,1>& dydk)
{
  double ak = Ak(k);
  double ck = Ck(k);
  double fk = (2.*k*k - 1.)*ak  - 2.*k*k*ck;

  if (k > 1.)
    {
      double ek = Ek(k);

      dydk(0) = 1./k/k/ek;
      dydk(1) = fk*fk /ak;
      dydk(2) = (1. - y(0)*k*ak*ek/ck) * (1./ck - ck/ak/ek);
      dydk(3) = y(0) * (1. - ck*ck/ak/ek) * k;
    }
  else
    {
      dydk(0) = 0.;
      dydk(1) = fk*fk /ak;
      dydk(2) = 0.;
      dydk(3) = 0.;
    }
}
 
// ####
// A(k)
// ####
double Integrals::Ak (double k)
{
  if (k > 1.)
    {
      double invk = 1. /k;
      double Kk   = gsl_sf_ellint_Kcomp (invk, GSL_PREC_DOUBLE);
  
      return (2./M_PI) * Kk;
    }
  else
    {
      return (2./M_PI) * k * gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
    }
}

// ####
// C(k)
// ####
double Integrals::Ck (double k)
{
  if (k > 1.)
    {
      double invk = 1. /k;
      double Ek   = gsl_sf_ellint_Ecomp (invk, GSL_PREC_DOUBLE);
  
      return (2./M_PI) * Ek;
    }
  else
    {
      return (2./M_PI) * (gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE) + (k*k-1.) * gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE)) /k;
    }
}

// ####
// E(k)
// ####
double Integrals::Ek (double k)
{
  double invk = 1. /k;
  double Kk   = gsl_sf_ellint_Kcomp (invk, GSL_PREC_DOUBLE);
  double Ek   = gsl_sf_ellint_Ecomp (invk, GSL_PREC_DOUBLE);
   
  return (2./3./M_PI) * (2. * (2. - 1./k/k) * Ek - (1. - 1./k/k) * Kk);
}

// ######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive fourth-order Runge-Kutta scheme
//
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step-length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... step-length cannot change by more than this factor from
//                  step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step-length
//     h_max   ... maximum allowable step-length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step whilst attempting to maintain 
//  constant truncation error per step of acc:
//
//    flag = 0 ... error is absolute
//    flag = 1 ... error is relative
//    flag = 2 ... error is mixed
//
//  If step-length falls below h_min then routine aborts
// ######################################################################
void Integrals::RK4Adaptive (double& x, Array<double,1>& y, double& h, 
			     double& t_err, double acc, double S, int& rept,
			     int maxrept, double h_min, double h_max, int flag, 
			     int diag, FILE* file)
{
  int neqns = y.extent(0);
  Array<double,1> y0(neqns), y1(neqns);
  double hin = h;

  // Save initial data
  double x0 = x;
  y0 = y;

  // Take full step 
  RK4Fixed (x, y, h);

  // Save data
  y1 = y;

  // Restore initial data 
  x = x0;
  y = y0;

  // Take two half-steps 
  RK4Fixed (x, y, h/2.);
  RK4Fixed (x, y, h/2.);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (y(i) - y1(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs ((y(i) - y1(i)) / y(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = fabs ((y(i) - y1(i)) / y(i));
          err2  =  fabs (y(i) - y1(i));
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15) t_err = 1.e-15;

  // Calculate new step-length 
  double h_est = h * pow (fabs (acc / t_err), 0.2);

  // Prevent step-length from changing by more than factor S
  if (h_est / h > S)
    h *= S;
  else if (h_est / h < 1. / S)
    h /= S;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Abort if step-length falls below h_min
  if (fabs(h) < h_min)
    { 
      //printf ("Integrals::RK4Adpative: Warning - |h| < hmin at x = %11.4e\n", x);
      //exit (1);
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %11.4e hin = %11.4e err = %11.4e acc = %11.4e hest = %11.4e hout = %11.4e count = %3d\n", 
	     x, hin, t_err, acc, h_est, h, count);

  // If truncation error acceptable take step 
  if ((t_err <= acc) || (count >= maxrept))
    {  
      rept  = count;
      count = 0;
    }
  // If truncation error unacceptable repeat step 
  else 
    {
      count++;
      x = x0;
      y = y0;
      RK4Adaptive (x, y, h, t_err, acc, S, rept, 
		   maxrept, h_min, h_max, flag, diag, file);
    }
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step-length fourth-order Runge-Kutta scheme.
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step-length
// #####################################################################
void Integrals::RK4Fixed (double& x, Array<double,1>& y, double h)
{
  int neqns = y.extent(0);
  Array<double,1> dydx(neqns), k1(neqns), k2(neqns), k3(neqns);
  Array<double,1> k4(neqns), f(neqns);

  // Integralsth intermediate step 
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1(i) = h * dydx(i);
      f(i)  = y(i) + k1(i) / 2.;
    }

  // First intermediate step 
  Rhs (x + h / 2., f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2(i) = h * dydx(i);
      f(i)  = y(i) + k2(i) / 2.;
    }

  // Second intermediate step 
  Rhs (x + h / 2., f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3(i) = h * dydx(i);
      f(i)  = y(i) + k3(i);
    }

  // Third intermediate step 
  Rhs (x + h, f, dydx);
  for (int i = 0; i < neqns; i++)
    k4(i) = h * dydx(i);

  // Actual step 
  for (int i = 0; i < neqns; i++)
    y(i) += k1(i) / 6. + k2(i) / 3. + k3(i) / 3. + k4(i) / 6.;
  x += h;
}

// #####################
// Function to open file
// #####################
FILE* Integrals::OpenFile (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Integrals::OpenFile: Error opening data-file\n");
      exit (1);
    }
  return file;
}
