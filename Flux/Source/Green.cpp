// Green.cpp

// PROGRAM ORGANIZATION:
//
// double Flux:: GreenPlasmaCos    (int i, int j, int ip, int jp)
// double Flux:: GreenPlasmaSin    (int i, int j, int ip, int jp)
// double Flux:: GreenPlasmaCosNeo (int i, int j, int ip, int jp)
// double Flux:: GreenPlasmaSinNeo (int i, int j, int ip, int jp)

#include "Flux.h"

// ########################################
// Functions to calculate Green's functions
// ########################################
double Flux::GreenPlasmaCos (int i, int j, int ip, int jp)
{
  double m      = double (mres[i]);
  double mp     = double (mres[ip]);
  double xheta  = th[j];
  double xhetap = th[jp];
  
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
  fun          *= cos (m*xheta - mp*xhetap);

  return fun;
}

double Flux::GreenPlasmaSin (int i, int j, int ip, int jp)
{
  double m      = double (mres[i]);
  double mp     = double (mres[ip]);
  double xheta  = th[j];
  double xhetap = th[jp];
  
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
  fun          *= - sin (m*xheta - mp*xhetap);

  return fun;
}

double Flux::GreenPlasmaCosNeo (int i, int j, int ip, int jp)
{
  double m      = double (mres[i]);
  double mp     = double (mres[ip]);
  double xheta  = gsl_matrix_get (theta, i,  j);
  double xhetap = gsl_matrix_get (theta, ip, jp);
  
  double R      = gsl_matrix_get (Rnc, i,  j);
  double Z      = gsl_matrix_get (Znc, i,  j);
  double Rp     = gsl_matrix_get (Rnc, ip, jp);
  double Zp     = gsl_matrix_get (Znc, ip, jp);

  double trans  = gsl_matrix_get (factor, i, j) * gsl_matrix_get (factor, ip, jp);

  double fac    = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp) + ETA);
  double eta    = atanh (fac);
  double ceta   = cosh (eta);

  double fun    = cos (double (NTOR) * M_PI) * M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  fun          *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  fun          *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  fun          *= cos (m*xheta - mp*xhetap);
  fun          *= trans;

  return fun;
}

double Flux::GreenPlasmaSinNeo (int i, int j, int ip, int jp)
{
  double m      = double (mres[i]);
  double mp     = double (mres[ip]);
  double xheta  = gsl_matrix_get (theta, i,  j);
  double xhetap = gsl_matrix_get (theta, ip, jp);
  
  double R      = gsl_matrix_get (Rnc, i,  j);
  double Z      = gsl_matrix_get (Znc, i,  j);
  double Rp     = gsl_matrix_get (Rnc, ip, jp);
  double Zp     = gsl_matrix_get (Znc, ip, jp);

  double trans  = gsl_matrix_get (factor, i, j) * gsl_matrix_get (factor, ip, jp);

  double fac    = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp)*(Z - Zp) + ETA);
  double eta    = atanh (fac);
  double ceta   = cosh (eta);

  double fun    = cos (double (NTOR) * M_PI) *  M_PI*M_PI * R*Rp /2. /gsl_sf_gamma (0.5) /gsl_sf_gamma (double (NTOR) + 0.5);
  fun          *= sqrt (ceta / (R*R + Rp*Rp + (Z - Zp)*(Z - Zp)));
  fun          *= (double (NTOR) - 0.5) * ToroidalP (NTOR-1, 0, ceta) + ToroidalP (NTOR+1, 0, ceta) /(double (NTOR) + 0.5);
  fun          *= - sin (m*xheta - mp*xhetap);
  fun          *= trans;

  return fun;
}

