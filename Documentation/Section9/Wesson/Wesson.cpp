// Wesson.cpp

#include "Wesson.h"

// ###########
// Constructor
// ###########
Wesson::Wesson ()
{
  // --------------------------------------
  // Set default values of class parameters
  // --------------------------------------
  mpol = 2.;
  ntor = 1.;
  rw   = 1.1;
  q0   = 2.;
  qa   = 3.;
  
  eps   = 1.e-6;
  delta = 1.e-8;
  h0    = 1.e-6;
  acc   = 1.e-13;
  hmin  = 1.e-10;
}

// #########################
// Function to solve problem
// #########################
void Wesson::Solve ()
{
  double two  = 2.;
  double one  = 1.;
  double zero = 0.;

  q0   = 0.8;
  mpol = 2.;
  ntor = 1.;
  rw   = 1.2;
  
  int    N       = 1000;
  double qastart = 2.0001;
  double qaend   = 20.;

  FILE* file = fopen ("m2n1r12.out", "w");
  fprintf (file, "%11.4e %11.4e %11.4e %11.4e\n", two, one, zero, zero);
  for (int i = 0; i < N; i++)
    {
      qa = qastart + (double (i) / double (N)) * (qaend - qastart);

      double Deltas = GetDelta ();

      double rs     = Findrs (mpol /ntor);
      double alphas = Getalpha (rs);
      double betas  = Getbeta (rs);
      double Wsat;
      if (Deltas > 0.)
	Wsat = rs * Deltas / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas);
      else
	Wsat = 0.;
      
      printf ("qa = %11.4e  rs = %11.4e  Delta = %11.4e  Wsat = %11.4e\n", qa, rs,  Deltas, Wsat);
      fprintf (file, "%11.4e %11.4e %11.4e %11.4e\n", qa, rs, Deltas, Wsat);
    }
  fclose (file);

  rw = 1.1;
  file = fopen ("m2n1r11.out", "w");
  fprintf (file, "%11.4e %11.4e %11.4e %11.4e\n", two, one, zero, zero);
  for (int i = 0; i < N; i++)
    {
      qa = qastart + (double (i) / double (N)) * (qaend - qastart);

      double Deltas = GetDelta ();

      double rs     = Findrs (mpol /ntor);
      double alphas = Getalpha (rs);
      double betas  = Getbeta (rs);
      double Wsat;
      if (Deltas > 0.)
	Wsat = rs * Deltas / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas);
      else
	Wsat = 0.;
      
      printf ("qa = %11.4e  rs = %11.4e  Delta = %11.4e  Wsat = %11.4e\n", qa, rs,  Deltas, Wsat);
      fprintf (file, "%11.4e %11.4e %11.4e %11.4e\n", qa, rs, Deltas, Wsat);
    }
  fclose (file);

  rw = 1.0001;
  file = fopen ("m2n1r10.out", "w");
  fprintf (file, "%11.4e %11.4e %11.4e %11.4e\n", two, one, zero, zero);
  for (int i = 0; i < N; i++)
    {
      qa = qastart + (double (i) / double (N)) * (qaend - qastart);

      double Deltas = GetDelta ();

      double rs     = Findrs (mpol /ntor);
      double alphas = Getalpha (rs);
      double betas  = Getbeta (rs);
      double Wsat;
      if (Deltas > 0.)
	Wsat = rs * Deltas / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas);
      else
	Wsat = 0.;
      
      printf ("qa = %11.4e  rs = %11.4e  Delta = %11.4e  Wsat = %11.4e\n", qa, rs,  Deltas, Wsat);
      fprintf (file, "%11.4e %11.4e %11.4e %11.4e\n", qa, rs, Deltas, Wsat);
    }
  fclose (file);
}

// #############################################
// Function to calculate tearing stability index
// #############################################
double Wesson::GetDelta ()
{
  // Find radius of resonant surface
  double qs = mpol /ntor;
  double rs = Findrs (qs);

  // Launch solution from magnetic axis and integrate to resonant surface
  double r, h, t_err;
  int rept; count = 0;
  Array<double,1> y(2);

  double nu    = Getnu ();
  double kappa = - nu / (1. - q0 /qs) /(mpol + 1.);

  r    = eps;
  y(0) = pow (r, mpol) + kappa * pow (r, mpol + 2);
  y(1) = mpol * pow (r, mpol - 1.) + kappa * (mpol + 2.) * pow (r, mpol + 1.);
  h    = h0;

  do
    {
      RK4Adaptive (r, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-2, 2, 0, NULL);
    }
  while (r < rs - delta);
  RK4Fixed (r, y, rs - delta - r);

  double alphas = Getalpha (r);
  double Deltam = r * y(1) /y(0) - alphas * (1. + log(delta));

  //printf ("rm = %11.4e  rs - rm = %11.4e  Deltam = %11.4e\n", r, rs - r, Deltam);

  // Launch solution from plasma boudnary and integrate to resonant surface

  double fac = pow (1./rw, 2.*mpol);
  r    = 1.;
  y(0) = 1.; 
  y(1) = - mpol * (1. + fac) /(1. - fac);
  h    = - h0;

  do
    {
      RK4Adaptive (r, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-2, 2, 0, NULL);
    }
  while (r > rs + delta);
  RK4Fixed (r, y, rs + delta - r);

  alphas = Getalpha (r);
  double Deltap = r * y(1) /y(0) - alphas * (1. + log(delta));

  //printf ("rm = %11.4e  rs - rm = %11.4e  Deltap = %11.4e\n", r, rs - r, Deltap);

  return Deltap - Deltam;
}

// ########################################
// Function to find resonant surface radius
// ########################################
double Wesson::Findrs (double qs)
{
  if (qs <= q0 || qs >= qa)
    {
      printf ("Wesson:Findrs - Error: no resonant surface in plasma\n");
      exit (1); 
    }

  int N = 1000;
  for (int i = 0; i < N; i++)
    {
      double r1 = double (i) /double (N);
      double r2 = double (i+1) /double (N);
      double q1 = Getq (r1);
      double q2 = Getq (r2);
      //printf ("i = %3d  r1 = %11.4e  q1 = %11.4e  r2 = %11.4e  q2 = %11.4e\n", i, r1, q1, r2, q2);
      
      if ((q1 - qs) * (q2 - qs) < 0.)
	{
	  if (i == 0)
	    {
	      double nu = Getnu ();
	      return sqrt ((qs/q0 - 1.) * (nu/2.));
	    }
	  else
	    {
	      double rn = r1, qn;
	      int count = 0;
	      do
		{
		  qn = Getq (rn);
		  double sn = Gets (rn);
		  
		  //printf ("rn = %11.4e  qn = %11.4e  residual = %11.4e\n", rn, qn, fabs (qn - qs));
		  
		  rn += rn * (qs /qn - 1.) /sn;
		  count++;
		}
	      while (fabs (qn - qs) > 1.e-15 && count < 100);
	      
	      if (isnan (rn))
		{
		  printf ("Wesson:Findrs - Error: rn is NaN\n");
		  exit (1);
		}
	      else
		return rn;
	    }
	}
    }
}

// ########################################
// Function to return value of parameter nu
// ########################################
double Wesson::Getnu ()
{
  return qa /q0 - 1.;
}

// ##########################################
// Function to return value of plasma current
// ##########################################
double Wesson::GetJ (double r)
{
  double nu = Getnu ();
  double r2 = r*r;

  return 2. * pow (1. - r2, nu) / q0;
}

// ##############################################################
// Function to return value of first derivative of plasma current
// ##############################################################
double Wesson::GetJp (double r)
{
  double nu = Getnu ();
  double r2 = r*r;

  return - 4. * nu * r * pow (1. - r2, nu - 1.) /q0;
}

// ###############################################################
// Function to return value of second derivative of plasma current
// ###############################################################
double Wesson::GetJpp (double r)
{
  double nu = Getnu ();
  double r2 = r*r;

  return - 4. * nu * pow (1. - r2, nu - 2.) * (1. - (2.*nu - 1.) * r2) /q0;
}

// #########################################
// Function to return value of safety-factor
// #########################################
double Wesson::Getq (double r)
{
  double nu = Getnu ();
  double r2 = r*r;

  if (r < eps)
    return q0 * (1. + 0.5 * nu * r2);
  else
    return qa * r2 /(1. - pow (1. - r2, nu + 1.));
}

// ##########################################
// Function to return value of magnetic shear
// ##########################################
double Wesson::Gets (double r)
{
  double J = GetJ (r);
  double q = Getq (r);

  return 2. - q * J; 
}

// ###########################################
// Function to return value of parameter alpha
// ###########################################
double Wesson::Getalpha (double r)
{
  double Jp = GetJp (r);
  double q  = Getq  (r);
  double s  = Gets  (r);

  return - q * r * Jp /s; 
}

// ##########################################
// Function to return value of parameter beta
// ##########################################
double Wesson::Getbeta (double r)
{
  double Jpp = GetJpp (r);
  double q   = Getq   (r);
  double s   = Gets   (r);

  return - q * r*r * Jpp /s; 
}

// ##############################################################
// Function to evaluate right-hand side of differential equations
// ##############################################################
void Wesson::Rhs (double r, Array<double,1>& y, Array<double,1>& dydr)
{
  dydr(0) = y(1);
  dydr(1) = - y(1) /r + mpol*mpol * y(0) /r/r + GetJp (r) * y(0) /r /(1./Getq (r) - ntor /mpol);
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
void Wesson::RK4Adaptive (double& x, Array<double,1>& y, double& h, 
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
      //printf ("Wesson::RK4Adpative: Warning - |h| < hmin at x = %11.4e\n", x);
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
void Wesson::RK4Fixed (double& x, Array<double,1>& y, double h)
{
  int neqns = y.extent(0);
  Array<double,1> dydx(neqns), k1(neqns), k2(neqns), k3(neqns);
  Array<double,1> k4(neqns), f(neqns);

  // Wessonth intermediate step 
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
FILE* Wesson::OpenFile (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Wesson::OpenFile: Error opening data-file\n");
      exit (1);
    }
  return file;
}
