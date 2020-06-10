// Green.cpp

// PROGRAM ORGANIZATION:
//
// double Flux:: GreenPlasmaCos   (int i, int j, int ip, int jp)
// double Flux:: GreenPlasmaSin   (int i, int j, int ip, int jp)
// double Flux:: GreenInboardCos  (int i, int j)
// double Flux:: GreenInboardSin  (int i, int j)
// double Flux:: GreenOutboardCos (int i, int j)
// double Flux:: GreenOutboardSin (int i, int j)

#include "Flux.h"

// ########################################
// Functions to calculate Green's functions
// ########################################
double Flux::GreenPlasmaCos (int i, int j, int ip, int jp)
{
  double m      = double (mres [i]);
  double mp     = double (mres [ip]);
  double theta  = th [j];
  double thetap = th [jp];
  
  double R      = gsl_matrix_get (Rst, i,  j);
  double Z      = gsl_matrix_get (Zst, i,  j);
  double Rp     = gsl_matrix_get (Rst, ip, jp);
  double Zp     = gsl_matrix_get (Zst, ip, jp);

  double fac    = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp) + ETA);
  double eta    = atanh (fac);
  double ceta   = cosh (eta);

  double fun    = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  fun          *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  fun          *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  fun          *= cos (m*theta - mp*thetap);

  return fun;
}

double Flux::GreenPlasmaSin (int i, int j, int ip, int jp)
{
  double m      = double (mres [i]);
  double mp     = double (mres [ip]);
  double theta  = th [j];
  double thetap = th [jp];
  
  double R      = gsl_matrix_get (Rst, i,  j);
  double Z      = gsl_matrix_get (Zst, i,  j);
  double Rp     = gsl_matrix_get (Rst, ip, jp);
  double Zp     = gsl_matrix_get (Zst, ip, jp);

  double fac    = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp) + ETA);
  double eta    = atanh (fac);
  double ceta   = cosh (eta);

  double fun    = cos (double (NTOR) * M_PI) *  M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  fun          *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  fun          *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  fun          *= - sin (m*theta - mp*thetap);

  return fun;
}

double Flux::GreenInboardCos (int i, int j)
{
  double m     = double (mres [i]);
  double theta = th [j];

  double R     = gsl_matrix_get (Rst, i, j);
  double Z     = gsl_matrix_get (Zst, i, j);

  double Rp    = RIN /R0 + DR;
  double Zp    = Zaxis;

  double fac   = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  double eta   = atanh (fac);
  double ceta  = cosh (eta);

  double funp  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funp        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funp        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funp        *= cos (m*theta);

  Rp           = RIN /R0 - DR;
  Zp           = Zaxis;

  fac          = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  eta          = atanh (fac);
  ceta         = cosh (eta);

  double funm  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funm        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funm        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funm        *= cos (m*theta);

  return (funm - funp) /2./DR;
}

double Flux::GreenInboardSin (int i, int j)
{
  double m     = double (mres [i]);
  double theta = th [j];

  double R     = gsl_matrix_get (Rst, i, j);
  double Z     = gsl_matrix_get (Zst, i, j);

  double Rp    = RIN /R0 + DR;
  double Zp    = Zaxis;

  double fac   = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  double eta   = atanh (fac);
  double ceta  = cosh (eta);

  double funp  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funp        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funp        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funp        *= - sin (m*theta);

  Rp           = RIN /R0 - DR;
  Zp           = Zaxis;

  fac          = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  eta          = atanh (fac);
  ceta         = cosh (eta);

  double funm  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funm        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funm        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funm        *= - sin (m*theta);

  return (funm - funp) /2./DR;
}

double Flux::GreenOutboardCos (int i, int j)
{
  double m     = double (mres [i]);
  double theta = th [j];

  double R     = gsl_matrix_get (Rst, i, j);
  double Z     = gsl_matrix_get (Zst, i, j);

  double Rp    = ROUT /R0 + DR;
  double Zp    = Zaxis;

  double fac   = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  double eta   = atanh (fac);
  double ceta  = cosh (eta);

  double funp  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funp        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funp        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funp        *= cos (m*theta);

  Rp           = ROUT /R0 - DR;
  Zp           = Zaxis;

  fac          = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  eta          = atanh (fac);
  ceta         = cosh (eta);

  double funm  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funm        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funm        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funm        *= cos (m*theta);

  return (funm - funp) /2./DR;
}

double Flux::GreenOutboardSin (int i, int j)
{
  double m     = double (mres [i]);
  double theta = th [j];

  double R     = gsl_matrix_get (Rst, i, j);
  double Z     = gsl_matrix_get (Zst, i, j);

  double Rp    = ROUT /R0 + DR;
  double Zp    = Zaxis;

  double fac   = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  double eta   = atanh (fac);
  double ceta  = cosh (eta);

  double funp  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funp        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funp        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funp        *= - sin (m*theta);

  Rp           = ROUT /R0 - DR;
  Zp           = Zaxis;

  fac          = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp));
  eta          = atanh (fac);
  ceta         = cosh (eta);

  double funm  = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  funm        *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  funm        *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  funm        *= - sin (m*theta);

  return (funm - funp) /2./DR;
}

