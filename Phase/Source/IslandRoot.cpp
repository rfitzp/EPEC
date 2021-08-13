// IslandRoot.cpp

// #####################
// PROGRAM ORGANIZATION:
// #####################

// double Phase::GetIslandWidth       (int j)
// double Phase::GetVacuumIslandWidth (int j)
// double Phase::GetVacuumIslandWidth (int j, double chi)
// void   Phase::GetIslandLimits      (int j, double Psi, double& Xminus, double& Xplus)
// double Phase::GetIslandRoot        (double c)

#include "Phase.h"

// ############################################################
// Function to find full width in PsiN of magnetic island chain
// ############################################################
double Phase::GetIslandWidth (int j)
{
  double Xminus, Xplus;

  GetIslandLimits (j, Psik (j), Xminus, Xplus);

  return Xplus - Xminus;
}

// ###################################################################
// Function to find full width in PsiN of vacuum magnetic island chain
// ###################################################################
double Phase::GetVacuumIslandWidth (int j)
{
  double Xminus, Xplus;

  GetIslandLimits (j, chi (j), Xminus, Xplus);

  return Xplus - Xminus;
}

// ###################################################################
// Function to find full width in PsiN of vacuum magnetic island chain
// ###################################################################
double Phase::GetVacuumIslandWidth (int j, double chi)
{
  double Xminus, Xplus;

  GetIslandLimits (j, chi, Xminus, Xplus);

  return Xplus - Xminus;
}

// ########################################################
// Function to find limits in PsiN of magnetic island chain 
// ########################################################
void Phase::GetIslandLimits (int j, double Psi, double& Xminus, double& Xplus)
{
  double cplus  = 2. * A1 (j) * fabs (Psi) /Deltakp (j)/Deltakp (j);
  double cminus = 2. * A1 (j) * fabs (Psi) /Deltakm (j)/Deltakm (j);
  
  Xplus  = + Deltakp (j) * GetIslandRoot (cplus);
  Xminus = - Deltakm (j) * GetIslandRoot (cminus);
}

// ############################################################
// Function to solve - x - log (1 - x) = c via Newton iteration
// ############################################################
double Phase::GetIslandRoot (double c)
{
  int    NITER = 10;
  double XMAX  = 0.999;
  
  double x = 0.5;
  
  for (int i = 0; i < NITER; i++)
    {
      double f = - x - log (1. - x);
      double fp = x /(1. - x);

      x += (c - f) /fp;

      if (x > 1.)
	x = XMAX;
    }

  return x;
}
