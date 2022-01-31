// Integrate.cpp

// PROGRAM ORGANIZATION:
//
// void Flux:: CalcQGP                        ()
// void Flux:: CheckQP                        ()
// void Flux:: CalcrP                         ()
// void FLux:: CalcGGJ                        ()
// void Flux:: CalcStraightAngle              ()
// void Flux:: CalcGamma                      ()
// void Flux:: CalcNeoclassicalAngle          ()
// void Flux:: CalcTearingSolutionNoWall      (int i)
// void Flux:: CalcTearingSolutionPerfectWall (int i)
// int  Flux:: Rhs1                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs2                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs3                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs4                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs5                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs6                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs7                           (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs8                           (double r, const double y[], double dydr[], void*)

#include "Flux.h"

// #######################################
// Function to calculate q(P)/g(P) profile
// #######################################
void Flux::CalcQGP ()
{
  gsl_odeiv2_system           sys1 = {pRhs1, NULL, 4, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys1, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [4]; 
  double                      r;
  
  for (int j = 1; j < NPSI; j++)
    {
      r    = 0.;
      y[0] = RP[j];
      y[1] = ZPTS[jc];
      y[2] = 0.;
      y[3] = 0.;
  
      while (r < M_PI)
	{
	  int status = gsl_odeiv2_driver_apply (d, &r, M_PI, y);

	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}
      RP1[j] = y[0];

      while (r < 2.*M_PI)
	{
	  int status = gsl_odeiv2_driver_apply (d, &r, 2.*M_PI, y);

	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}
      
      QGP[j] = y[2];
      QP [j] = QGP[j]*GP[j];
      J0 [j] = y[3];

      if (j%10 == 0 || j > NPSI-10)
	{
	  printf ("j = %4d  PsiN = %11.4e  q = %11.4e\n", j, 1.-P[j], QP[j]);
	  fflush (stdout);
	}
    }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ################################################
// Function to confirm q value on rational surfaces
// ################################################
void Flux::CheckQP ()
{
  gsl_odeiv2_system           sys1 = {pRhs1, NULL, 4, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys1, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [4]; 
  double                      r;
  
  for (int j = 0; j < nres; j++)
    {
      r    = 0.;
      y[0] = Rres[j];
      y[1] = ZPTS[jc];
      y[2] = 0.;
      y[3] = 0.;
  
      while (r < 2.*M_PI)
	{
	  int status = gsl_odeiv2_driver_apply (d, &r, 2.*M_PI, y); 
	  
	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}

      double qgp = y[2];
      double qp  = qgp * gres[j];
      printf ("mpol = %3d  qres = %11.4e  q = %11.4e  residual = %11.4e\n",
	      mres[j], double (mres[j]) /double(NTOR), qp, 1. - double (mres[j]) /double(NTOR)/qp);
      fflush (stdout);
    }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ##################################
// Function to calculate r(P) profile
// ##################################
void Flux::CalcrP ()
{
  gsl_odeiv2_system           sys2 = {pRhs2, NULL, 1, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys2, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [1]; 
  double                      r;

  // Calculate r[Psi]
  for (int j = 1; j < NPSI; j++)
    {
      r    = P[j];
      y[0] = 0.;
      
      while (r < P[0])
	{
	  int status = gsl_odeiv2_driver_apply (d, &r, P[0], y);
	  
	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}

      rP[j] = sqrt (y[0]);
    }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ##############################
// Function to calculate GGJ data
// ##############################
void Flux::CalcGGJ ()
{
  gsl_odeiv2_system           sys6 = {pRhs6, NULL, 8, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys6, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [8]; 
  double                      r;
  
  for (int j = 0; j < nres; j++)
    {
      r    = 0.;
      y[0] = Rres[j];
      y[1] = ZPTS[jc];
      y[2] = 0.;
      y[3] = 0.;
      y[4] = 0.;
      y[5] = 0.;
      y[6] = 0.;
      y[7] = 0.;
      qgp  = gres[j];
  
      while (r < 2.*M_PI)
	{
	  int status =  gsl_odeiv2_driver_apply (d, &r, 2.*M_PI, y); 

	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}
      
      J1[j] = y[2];
      J2[j] = y[3];
      J3[j] = y[4];
      J4[j] = y[5];
      J5[j] = y[6];
      J6[j] = y[7];
  }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// #########################################
// Function to calculate straight angle data
// #########################################
void Flux::CalcStraightAngle ()
{
  gsl_odeiv2_system           sys3 = {pRhs3, NULL, 2, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys3, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [2]; 
  double                      r, psir, psiz, rt, zt;
  
  for (int j = 0; j < nres; j++)
    {
      r    = 0.;
      y[0] = Rres[j];
      y[1] = ZPTS[jc];
      qgp  = fabs (Psic) * qres[j] /gres[j];
	      
      gsl_matrix_set (Rst, j, 0, y[0]);
      gsl_matrix_set (Zst, j, 0, y[1]);

      for (int k = 1; k < NTHETA; k++)
	{
	  while (r < th[k])
	    {
	      int status = gsl_odeiv2_driver_apply (d, &r, th[k], y); 
	      
	      if (status != GSL_SUCCESS)
		{
		  printf ("Profile: status != GSL_SUCCESS\n");
		  exit (1);
		}
	    }
	  gsl_matrix_set (Rst, j, k, y[0]);
	  gsl_matrix_set (Zst, j, k, y[1]);
	}
     }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// #################################################
// Function to calculate gammas on rational surfaces
// #################################################
void Flux::CalcGamma ()
{
  gsl_odeiv2_system           sys4 = {pRhs4, NULL, 3, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys4, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [3]; 
  double                      r;
  
  for (int j = 0; j < nres; j++)
    {
      r    = 0.;
      y[0] = Rres[j];
      y[1] = ZPTS[jc];
      y[2] = 0.;
      qgp  = gres[j];

      while (r < 2.*M_PI)
	{
	  int status = gsl_odeiv2_driver_apply (d, &r, 2.*M_PI, y); 

	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}

      gmres[j] = 1./y[2];
    }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// #############################################
// Function to calculate neoclassical angle data
// #############################################
void Flux::CalcNeoclassicalAngle ()
{
  gsl_odeiv2_system           sys5 = {pRhs5, NULL, 3, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys5, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [3]; 
  double                      r, psir, psiz, rt, zt;
  
  for (int j = 0; j < nres; j++)
    {
      r    = 0.;
      y[0] = Rres[j];
      y[1] = ZPTS[jc];
      y[2] = 0.;
      qgp  = gres[j];
      qgp1 = fabs (Psic) /gmres[j];
      qgp2 = gres[j] /gmres[j] /qres[j];
	      
      gsl_matrix_set (Rnc,   j, 0, y[0]);
      gsl_matrix_set (Znc,   j, 0, y[1]);
      gsl_matrix_set (theta, j, 0, y[2]);

      for (int k = 1; k < NTHETA; k++)
	{
	  while (r < Th[k])
	    {
	      int status =  gsl_odeiv2_driver_apply (d, &r, Th[k], y); 
	      
	      if (status != GSL_SUCCESS)
		{
		  printf ("Profile: status != GSL_SUCCESS\n");
		  exit (1);
		}
	    }
	  gsl_matrix_set (Rnc,   j, k, y[0]);
	  gsl_matrix_set (Znc,   j, k, y[1]);
	  gsl_matrix_set (theta, j, k, y[2]);
	}
     }
  
  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ###################################################################################
// Function to solve no wall cylindrical tearing mode problem for ith rational surface
// ###################################################################################
void Flux::CalcTearingSolutionNoWall (int i)
{
  gsl_odeiv2_system           sys7 = {pRhs7, NULL, 2, this};
  gsl_odeiv2_system           sys8 = {pRhs8, NULL, 2, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys7, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [2]; 
  double                      r, c, rin, rout;

  ires = i;

  double rs = rres[i];
  double m  = double (mres[i]);
  double a  = alpres[i];
  double b  = betres[i];
  double D  = DELTA;
  double lD = log (DELTA);
  
  r    = EPS;
  c    = Interpolate (NPSI, rP, JPrr, r, 0) /4. /(1. + m);
  y[0] = 1. + c * r*r;
  y[1] = 2. * c * r;
  rin  = rs - D;
  
  int status = gsl_odeiv2_driver_apply (d, &r, rin, y); 
  
  if (status != GSL_SUCCESS)
    {
      printf ("Profile: status != GSL_SUCCESS\n");
      exit (1);
    }

  double a1m = a * (1. + lD) - (m*m /rs/rs + b - 0.5 * a*a * (3. - 1./a/rs) + a*a * (1. - 1./a/rs) * (0.5 + lD)) * D;
  double b1m = 0.5 - 0.5 * a * (2. - 1./a/rs) * D;
  double a2m = 1. - a * D * lD;
  double b2m = - 0.5 * D;
  
  double psim  = pow (r, m) * y[0];
  double psipm = m * pow (r, m-1.) * y[0] + pow (r, m) * y[1];

  double lm  = psipm /psim;
  double SDm = (a2m * lm - a1m) /(b1m - b2m * lm);
  
  gsl_odeiv2_driver_free (d);

  d    = gsl_odeiv2_driver_alloc_y_new (&sys8, T, - H0, ACC/2., ACC/2.);
  r    = ra;
  y[0] = 1.;
  y[1] = 0.;
  rout = rs + D;

  status = gsl_odeiv2_driver_apply (d, &r, rout, y); 
  
  if (status != GSL_SUCCESS)
    {
      printf ("Profile: status != GSL_SUCCESS\n");
      exit (1);
    }

  double a1p = a * (1. + lD) + (m*m /rs/rs + b - 0.5 * a*a * (3. - 1./a/rs) + a*a * (1. - 1./a/rs) * (0.5 + lD)) * D;
  double b1p = 0.5 + 0.5 * a * (2. - 1./a/rs) * D;
  double a2p = 1. + a * D * lD;
  double b2p = + 0.5 * D;
  
  double psip  = pow (r, -m) * y[0];
  double psipp = - m * pow (r, -m-1.) * y[0] + pow (r, -m) * y[1];

  double lp  = psipp /psip;
  double SDp = (a2p * lp - a1p) /(b1p - b2p * lp);
  
  Delta_nw[i] = rs * 0.5 * (SDp - SDm);
  Sigma_nw[i] = rs * 0.5 * (SDp + SDm);

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ########################################################################################
// Function to solve perfect wall cylindrical tearing mode problem for ith rational surface
// ########################################################################################
void Flux::CalcTearingSolutionPerfectWall (int i)
{
  gsl_odeiv2_system           sys7 = {pRhs7, NULL, 2, this};
  gsl_odeiv2_system           sys8 = {pRhs8, NULL, 2, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys7, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [2]; 
  double                      r, c, rin, rout;

  ires = i;

  double rs = rres[i];
  double m  = double (mres[i]);
  double a  = alpres[i];
  double b  = betres[i];
  double D  = DELTA;
  double lD = log (DELTA);
  
  r    = EPS;
  c    = Interpolate (NPSI, rP, JPrr, r, 0) /4. /(1. + m);
  y[0] = 1. + c * r*r;
  y[1] = 2. * c * r;
  rin  = rs - D;
  
  int status = gsl_odeiv2_driver_apply (d, &r, rin, y); 
  
  if (status != GSL_SUCCESS)
    {
      printf ("Profile: status != GSL_SUCCESS\n");
      exit (1);
    }

  double a1m = a * (1. + lD) - (m*m /rs/rs + b - 0.5 * a*a * (3. - 1./a/rs) + a*a * (1. - 1./a/rs) * (0.5 + lD)) * D;
  double b1m = 0.5 - 0.5 * a * (2. - 1./a/rs) * D;
  double a2m = 1. - a * D * lD;
  double b2m = - 0.5 * D;
  
  double psim  = pow (r, m) * y[0];
  double psipm = m * pow (r, m-1.) * y[0] + pow (r, m) * y[1];

  double lm  = psipm /psim;
  double SDm = (a2m * lm - a1m) /(b1m - b2m * lm);
  
  gsl_odeiv2_driver_free (d);

  d    = gsl_odeiv2_driver_alloc_y_new (&sys8, T, - H0, ACC/2., ACC/2.);
  r    = ra;
  y[0] = 1.;
  y[1] = - (2.*m /ra) /(pow(RW, 2.*m) - 1.);
  rout = rs + D;

  status = gsl_odeiv2_driver_apply (d, &r, rout, y); 
  
  if (status != GSL_SUCCESS)
    {
      printf ("Profile: status != GSL_SUCCESS\n");
      exit (1);
    }

  double a1p = a * (1. + lD) + (m*m /rs/rs + b - 0.5 * a*a * (3. - 1./a/rs) + a*a * (1. - 1./a/rs) * (0.5 + lD)) * D;
  double b1p = 0.5 + 0.5 * a * (2. - 1./a/rs) * D;
  double a2p = 1. + a * D * lD;
  double b2p = + 0.5 * D;
  
  double psip  = pow (r, -m) * y[0];
  double psipp = - m * pow (r, -m-1.) * y[0] + pow (r, -m) * y[1];

  double lp  = psipp /psip;
  double SDp = (a2p * lp - a1p) /(b1p - b2p * lp);
  
  Delta_pw[i] = rs * 0.5 * (SDp - SDm);
  Sigma_pw[i] = rs * 0.5 * (SDp + SDm);
  Sigma_w [i] = 2.*m * pow(rs/ra/RW, m) /(1. - pow (1./RW, 2.*m)) /y[0];

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// #####################################################
// Function to evaluate right-hand sides of q/g equation
// #####################################################
int Flux::Rhs1 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z
  // y[2] - q/g
  // y[3] - J0
  
  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);
  double Rc   = Raxis - y[0];
  double Zc   = y[1]  - Zaxis;
  double fac  = (Rc*Rc + Zc*Zc) /(- Zc*PsiZ + Rc*PsiR);

  dydr[0] = - PsiZ * fac;
  dydr[1] = + PsiR * fac;
  dydr[2] = (1./2./M_PI /fabs(Psic)) * fac /y[0];
  dydr[3] = y[0] * (1./2./M_PI /fabs(Psic)) * fac;
    
  return GSL_SUCCESS;
}

// ###################################################
// Function to evaluate right-hand sides of r equation
// ###################################################
int Flux::Rhs2 (double r, const double y[], double dydr[], void*)
{
  // y[0] - (rP)^2    

  double qg = Interpolate (NPSI, PsiN, QGP, 1.-r, 0);

  dydr[0] = 2.*fabs(Psic) * qg;

  return GSL_SUCCESS;
}

// #######################################################
// Function to evaluate right-hand sides of theta equation
// #######################################################
int Flux::Rhs3 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z

  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);

  dydr[0] = - qgp * y[0] * PsiZ;
  dydr[1] = + qgp * y[0] * PsiR;

  return GSL_SUCCESS;
}

// #######################################################
// Function to evaluate right-hand sides of gamma equation
// #######################################################
int Flux::Rhs4 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z
  // y[2] - 1/gamma

  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);
  double Grad = PsiR*PsiR + PsiZ*PsiZ;
  double RB   = sqrt (qgp*qgp + Psic*Psic * Grad);
  double Rc   = Raxis - y[0];
  double Zc   = y[1]  - Zaxis;
  double fac  = (Rc*Rc + Zc*Zc) /(- Zc*PsiZ + Rc*PsiR);

  dydr[0] = - PsiZ * fac;
  dydr[1] = + PsiR * fac;
  dydr[2] = (1./2./M_PI /fabs(Psic)) * RB * fac;
    
  return GSL_SUCCESS;
}

// #######################################################
// Function to evaluate right-hand sides of Theta equation
// #######################################################
int Flux::Rhs5 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z
  // y[2] - theta

  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);
  double Grad = PsiR*PsiR + PsiZ*PsiZ;
  double RB   = sqrt (qgp*qgp + Psic*Psic * Grad);

  dydr[0] = - qgp1 * PsiZ /RB;
  dydr[1] = + qgp1 * PsiR /RB;
  dydr[2] =   qgp2 /RB /y[0];

  return GSL_SUCCESS;
}

// #####################################################
// Function to evaluate right-hand sides of GGJ equation
// #####################################################
int Flux::Rhs6 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z
  // y[2] - J1
  // y[3] - J2
  // y[4] - J3
  // y[5] - J4
  // y[6] - J5
  // y[7] - J6
  
  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);
  double Rc   = Raxis - y[0];
  double Zc   = y[1]  - Zaxis;
  double fac  = (Rc*Rc + Zc*Zc) /(- Zc*PsiZ + Rc*PsiR);
  double np2  = Psic*Psic * (PsiR*PsiR + PsiZ*PsiZ);
  double b2   = (qgp*qgp + np2) /y[0]/y[0];
  double inp2 = 1. /np2;
  double ib2  = 1. /b2;

  dydr[0] = - PsiZ * fac;
  dydr[1] = + PsiR * fac;
  dydr[2] = y[0]              * (1./2./M_PI /fabs(Psic)) * fac;
  dydr[3] = y[0] * b2         * (1./2./M_PI /fabs(Psic)) * fac;
  dydr[4] = y[0] * ib2        * (1./2./M_PI /fabs(Psic)) * fac;
  dydr[5] = y[0]       * inp2 * (1./2./M_PI /fabs(Psic)) * fac;
  dydr[6] = y[0] * b2  * inp2 * (1./2./M_PI /fabs(Psic)) * fac;
  dydr[7] = y[0] * ib2 * inp2 * (1./2./M_PI /fabs(Psic)) * fac;
     
  return GSL_SUCCESS;
}

// #################################################################################
// Function to evaluate right-hand sides of tearing equation inside rational surface
// #################################################################################
int Flux::Rhs7 (double r, const double y[], double dydr[], void*)
{
  // y[0] - f     
  // y[1] - f'

  double m = double (mres[ires]);
  double n = double (NTOR);

  double q   = Interpolate (NPSI, rP, QP,  r, 0);
  double Jp  = Interpolate (NPSI, rP, JPr, r, 0);
  double fac = Jp /r /(1./q - n/m);

  dydr[0] = y[1];
  dydr[1] = - (1. + 2. * m) * y[1] /r + fac * y[0];
      
  return GSL_SUCCESS;
}

// ##################################################################################
// Function to evaluate right-hand sides of tearing equation outside rational surface
// ##################################################################################
int Flux::Rhs8 (double r, const double y[], double dydr[], void*)
{
  // y[0] - g     
  // y[1] - g'

  double m = double (mres[ires]);
  double n = double (NTOR);

  double q   = Interpolate (NPSI, rP, QP,  r, 0);
  double Jp  = Interpolate (NPSI, rP, JPr, r, 0);
  double fac = Jp /r /(1./q - n/m);

  dydr[0] = y[1];
  dydr[1] = - (1. - 2. * m) * y[1] /r + fac * y[0];
     
  return GSL_SUCCESS;
}
