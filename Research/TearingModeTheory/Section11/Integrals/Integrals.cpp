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

  xmax  = 20.;
  eps   = 1.e-7;
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

  double x, h, t_err;
  int rept; count = 0;
  Array<double,1> y(4);

  x    = eps;
  y(0) = 0.;
  y(1) = 0.;
  y(2) = 0.;
  y(3) = 0.;
  h    = h0;

  do
    {
      RK4Adaptive (x, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-1, 2, 0, NULL);
    }
  while (x < xmax);

  printf ("\nI0 = %11.4e  I1 = %11.4e  I2 = %11.4e  I3 = %11.4e\n\n", y(0), y(2), y(2), y(3));

  double K00i = y(0);
  double K01i = y(1);
  double K11i = y(2);

  double K00e = 1. + y(1);
  double K01e = 1. + y(2);
  double K11e = 2. + y(3);

  double mu00i = K00i;
  double mu01i = (5./2.) * K00i - K01i;
  double mu11i = K11i - 5. * K01i + (25./4.) * K00i;

  double mu00e = K00e;
  double mu01e = (5./2.) * K00e - K01e;
  double mu11e = K11e - 5. * K01e + (25./4.) * K00e;

  printf ("mu00i = %11.4e  mu01i = %11.4e  mu11i = %11.4e\n",   mu00i, mu01i, mu11i);
  printf ("mu00e = %11.4e  mu01e = %11.4e  mu11e = %11.4e\n\n", mu00e, mu01e, mu11e);

  double fnc  = mu01i /mu00i;
  double fbs1 = - ((sqrt(2.) + 13./4.) * mu00e - (3./2.) * mu01e) /(1. + sqrt(2.));
  double fbs2 = fnc * ((sqrt(2.) + 13./4.) * mu00e - (3./2.) * mu01e) /(1. + sqrt(2.));
  double fbs3 = fbs1;
  double fbs4 = ((sqrt(2.) + 13./4.) * mu01e - (3./2.) * mu11e) /(1. + sqrt(2.));

  printf ("fnc = %11.4e  fbs1 = %11.4e  fbs2 = %11.4e  fbs3 = %11.4e  fbs4 = %11.4e\n\n", fnc, fbs1, fbs2, fbs3, fbs4);
}

// ##############################################################
// Function to evaluate right-hand side of differential equations
// ##############################################################
void Integrals::Rhs (double x, Array<double,1>& y, Array<double,1>& dydx)
{
  dydx(0) = F(x) /x;
  dydx(1) = F(x);
  dydx(2) = F(x) * x;
  dydx(3) = F(x) * x*x;
}
 
// ######
// psi(x)
// ######
double Integrals::psi (double x)
{
  return gsl_sf_erf (x) - psip (x);
}

// #######
// psip(x)
// #######
double Integrals::psip (double x)
{
  return (2. /sqrt(M_PI)) * x * exp (-x*x);
}

// ####
// F(x)
// ####
double Integrals::F (double x)
{
  return ((1. - 1./2./x) * psi (x) + psip (x)) * exp (-x);
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
