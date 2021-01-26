// Stage2.cpp

// PROGRAM ORGANIATION:
//
// void Flux:: Stage2                      ()
// void Flux:: Stage2SetSimpsonWeights     ()
// void Flux:: Stage2ReadData              ()
// void Flux:: Stage2CalcQ                 ()
// void Flux:: Stage2FindRational          ()
// void Flux:: Stage2CalcStraightAngle     ()
// void Flux:: Stage2CalcNeoclassicalAngle ()
// void Flux:: Stage2CalcNeoclassicalPara  ()
// void Flux:: Stage2CalcMatrices          ()

#include "Flux.h"

// ####################################################
// Function to input Stage1 data and output Stage2 data
// ####################################################
void Flux::Stage2 ()
{
  // ..............................
  // Set weights for Simpson's rule
  // ..............................
  Stage2SetSimpsonWeights ();
  
  // ................................
  // Read data for Stage2 calculation
  // ................................
  Stage2ReadData ();
  
  // ..........................
  // Calculate Stage2 q profile
  // ..........................
  Stage2CalcQ ();

  // ......................
  // Find rational surfaces
  // ......................
  Stage2FindRational ();
  
  // ..................................................
  // Calculate straight angle data on rational surfaces
  // ..................................................
  Stage2CalcStraightAngle ();

  // ......................................................
  // Calculate neoclassical angle data on rational surfaces
  // ......................................................
  Stage2CalcNeoclassicalAngle ();
 
  // ......................................................
  // Calculate neoclassical parameters at rational surfaces
  // ......................................................
  Stage2CalcNeoclassicalPara ();
  
  // ............................
  // Calculate stability matrices
  // ............................
  Stage2CalcMatrices ();

  // ............
  // Output fFile
  // ............
  FILE* file = OpenFilew ((char*) "Outputs/fFile");
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e\n",
	   R0, B0, ra * R0, q95, r95 /ra, qlim, rlim /ra, QP[0], QP[NPSI-1], NPSI, NTOR, nres, PSILIM);
  for (int j = 0; j < NPSI; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     1. - P[j], rP[j] /ra, - ra * Interpolate (NPSI, rP, P, rP[j], 1));
  for (int i = 0; i < nres; i++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres[i], rres[i]/ra, sres[i], gres[i], gmres[i], Ktres[i], Kares[i], fcres[i], ajj[i], PsiNres[i], dPsidr[i], Khres[i], A1res[i]);
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      fprintf (file, "%d %d %16.9e %16.9e\n", mres[i], mres[j],
	       GSL_REAL (gsl_matrix_complex_get (FF, i, j)), GSL_IMAG (gsl_matrix_complex_get (FF, i, j)));
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      fprintf (file, "%d %d %16.9e %16.9e\n", mres[i], mres[j],
	       GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
  for (int i = 0; i < nres; i++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e\n", mres[i],
	     GSL_REAL (gsl_vector_complex_get (EI, i)), GSL_IMAG (gsl_vector_complex_get (EI, i)),
	     GSL_REAL (gsl_vector_complex_get (EO, i)), GSL_IMAG (gsl_vector_complex_get (EO, i)));
  fclose (file);

  if (INTG > 0)
    {
      char* filename = new char[MAXFILENAMELENGTH];
      sprintf (filename, "Outputs/fFiles/f.%d", int (TIME));
      file = OpenFilew (filename);
      fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e\n",
	       R0, B0, ra * R0, q95, r95 /ra, qlim, rlim/ra, QP[0], QP[NPSI-1], NPSI, NTOR, nres, PSILIM);
      for (int j = 0; j < NPSI; j++)
	fprintf (file, "%16.9e %16.9e %16.9e\n",
		 1. - P[j], rP[j] /ra, - ra * Interpolate (NPSI, rP, P, rP[j], 1));
      for (int i = 0; i < nres; i++)
	fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
		 mres[i], rres[i]/ra, sres[i], gres[i], gmres[i], Ktres[i], Kares[i], fcres[i], ajj[i], PsiNres[i], dPsidr[i], Khres[i], A1res[i]);
      for (int i = 0; i < nres; i++)
	for (int j = 0; j < nres; j++)
	  fprintf (file, "%d %d %16.9e %16.9e\n", i, j,
		   GSL_REAL (gsl_matrix_complex_get (FF, i, j)), GSL_IMAG (gsl_matrix_complex_get (FF, i, j)));
      for (int i = 0; i < nres; i++)
	for (int j = 0; j < nres; j++)
	  fprintf (file, "%d %d %16.9e %16.9e\n", i, j,
		   GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
      for (int i = 0; i < nres; i++)
	fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e\n", i,
		 GSL_REAL (gsl_vector_complex_get (EI, i)), GSL_IMAG (gsl_vector_complex_get (EI, i)),
		 GSL_REAL (gsl_vector_complex_get (EO, i)), GSL_IMAG (gsl_vector_complex_get (EO, i)));
      fclose (file);

      sprintf (filename, "f.%d", int (TIME));
      file = OpenFilea ((char*) "Outputs/fFiles/Index");
      fprintf (file, "%s %19.6e\n", filename, TIME);
      fclose (file);

      file = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
      fprintf (file, "Wrote fFile %s\n", filename);
      fclose (file);
      
      delete[] filename;
    }
  
  // ........
  // Clean up
  // ........
  delete[] RPTS;  delete[] ZPTS;
  delete[] RBPTS; delete[] ZBPTS;
  delete[] RLPTS; delete[] ZLPTS;

  delete[] PSIN; delete[] G;   delete[] Pr;
  delete[] GGp;  delete[] Prp; delete[] Q;

  delete[] s; delete[] Rs;

  delete[] P;    delete[] RP;  delete[] rP; delete[] GP;
  delete[] QGP;  delete[] QP;  delete[] PP; delete[] GPP;
  delete[] PPP;  delete[] S;   delete[] QX; delete[] RP1;
  delete[] Bt;   delete[] Bt1; delete[] Bp; delete[] Bp1;

  delete[] PsiN; delete[] QPN; delete[] A1;  
  
  delete[] mres;    delete[] qres;  delete[] rres;  delete[] sres;
  delete[] gres;    delete[] Rres;  delete[] gmres; delete[] fcres;
  delete[] PsiNres; delete[] Rres1; delete[] A1res; 

  delete[] th; 
  gsl_matrix_free (Rst); gsl_matrix_free (Zst);

  delete[] Th;
  gsl_matrix_free (Rnc); gsl_matrix_free (Znc); 

  gsl_matrix_free (Bnc); gsl_matrix_free (Cnc);
 
  delete[] I1;    delete[] I2;    delete[] I3;    
  delete[] Ktres; delete[] Kares; delete[] ajj; delete[] dPsidr;
  delete[] Khres;

  gsl_matrix_free (I4); gsl_matrix_free (I5); gsl_matrix_free (I6);
  
  gsl_matrix_complex_free (FF); gsl_matrix_complex_free (EE);
  gsl_vector_complex_free (EI);  gsl_vector_complex_free (EO);
}

// #############################################
// Function to assign weights for Simpson's rule
// #############################################
void Flux::Stage2SetSimpsonWeights ()
{
  double h  = 2.*M_PI   /double (NTHETA-1);  // Step length for angular integrals
         hh = 0.9999999 /double (NTHETA-1);  // Step length for fraction of circulating particles integrals

  Weight1D.resize (NTHETA);
  weight1D.resize (NTHETA);
  Weight2D.resize (NTHETA, NTHETA);
  for (int j = 0; j < NTHETA; j++)
    {
    	Weight1D (j) = 0.;
	weight1D (j) = 0.;

	for (int k = 0; k < NTHETA; k++)
	  Weight2D (j, k) = 0.;
    }
  for (int j = 0; j < NTHETA-2; j += 2)
    {
      Weight1D (j)   +=      h /3.;
      Weight1D (j+1) += 4. * h /3.;
      Weight1D (j+2) +=      h /3.;

      weight1D (j)   +=      hh /3.;
      weight1D (j+1) += 4. * hh /3.;
      weight1D (j+2) +=      hh /3.;
      
      for (int k = 0; k < NTHETA-2; k += 2)
	{
	  Weight2D (j,   k)   +=       h*h /9.;
	  Weight2D (j+1, k)   +=  4. * h*h /9.;
	  Weight2D (j+2, k)   +=       h*h /9.;
	  Weight2D (j,   k+1) +=  4. * h*h /9.;
	  Weight2D (j+1, k+1) += 16. * h*h /9.;
	  Weight2D (j+2, k+1) +=  4. * h*h /9.;
	  Weight2D (j,   k+2) +=       h*h /9.;
	  Weight2D (j+1, k+2) +=  4. * h*h /9.;
	  Weight2D (j+2, k+2) +=       h*h /9.;
	}
    }
}

// #############################################
// Function to read data for Stage2 calculations
// #############################################
void Flux::Stage2ReadData ()
{
  printf ("Reading data from gFile:\n");
  
  // ..............
  // Read R0 and B0
  // ..............
  FILE* file = OpenFiler ((char*) "Outputs/Stage1/R0B0.txt");
  if (fscanf (file, "%lf %lf", &R0, &B0) != 2)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/R0B0.txt\n");
      exit (1);
    }
  fclose (file);

  // ................
  // Read array sizes
  // ................
  file = OpenFiler ((char*) "Outputs/Stage1/Points.txt");
  if (fscanf (file, "%d %d %d %d", &NRPTS, &NZPTS, &NBPTS, &NLPTS) != 4)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Points.txt\n");
      exit (1);
    }
  fclose (file);

  // ..............
  // Read R, Z grid
  // ..............
  RPTS = new double[NRPTS];  // R array
  ZPTS = new double[NZPTS];  // Z array

  file = OpenFiler ((char*) "Outputs/Stage1/R.txt");
  for (int i = 0; i < NRPTS; i++)
    if (fscanf (file, "%lf", &RPTS[i]) != 1)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/R.txt\n");
	exit (1);
      }
  fclose (file);
  file = OpenFiler ((char*) "Outputs/Stage1/Z.txt");
  for (int j = 0; j < NZPTS; j++)
    if (fscanf (file, "%lf", &ZPTS[j]) != 1)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Z.txt\n");
	exit (1);
      }
  fclose (file);

  dR  = RPTS[1] - RPTS[0];
  dZ  = ZPTS[1] - ZPTS[0];
  dR2 = dR*dR;
  dZ2 = dZ*dZ;
  dR3 = dR*dR2;
  dZ3 = dZ*dZ2;

  // .......................
  // Read magnetic axis data
  // .......................
  file = OpenFiler ((char*) "Outputs/Stage1/Axis.txt");
  if (fscanf (file, "%lf %lf %lf", &Raxis, &Zaxis) != 2)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Axis.txt\n");
      exit (1);
    }	
  fclose (file);

  // ..............................
  // Read boundary and limiter data
  // ..............................
  RBPTS = new double[NBPTS];  // R coordinates of boundary
  ZBPTS = new double[NBPTS];  // Z coordinates of boundary

  file = OpenFiler ((char*) "Outputs/Stage1/Boundary.txt");
  for (int i = 0; i < NBPTS; i++)
    if (fscanf (file, "%lf %lf", &RBPTS[i], &ZBPTS[i]) != 2)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Boundary.txt\n");
	exit (1);
      }
  fclose (file);

  RLPTS = new double[NLPTS];  // R coordinates of limiter
  ZLPTS = new double[NLPTS];  // Z coordinates of limiter

  file = OpenFiler ((char*) "Outputs/Stage1/Limiter.txt");
  for (int i = 0; i < NLPTS; i++)
    if (fscanf (file, "%lf %lf", &RLPTS[i], &ZLPTS[i]) != 2)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Limiter.txt\n");
	exit (1);
      }
  fclose (file);

  // .................
  // Find RIN and ROUT
  // .................
  for (int i = 0; i < NLPTS-1; i++)
    {
      if ((ZLPTS[i] - Zaxis) * (ZLPTS[i+1] - Zaxis) < 0.)
	{
	  if (RLPTS[i] < Raxis)
	    RIN  = (RLPTS[i] * (ZLPTS[i+1] - Zaxis) + RLPTS[i+1] * (Zaxis - ZLPTS[i])) /(ZLPTS[i+1] - ZLPTS[i]);
	  else if (RLPTS[i] > Raxis)
	    ROUT = (RLPTS[i] * (ZLPTS[i+1] - Zaxis) + RLPTS[i+1] * (Zaxis - ZLPTS[i])) /(ZLPTS[i+1] - ZLPTS[i]);
	}
    }
  if ((ZLPTS[NLPTS-1] - Zaxis) * (ZLPTS[0] - Zaxis) < 0.)
    {
      if (RLPTS[NLPTS-1] < Raxis)
	RIN  = (RLPTS[NLPTS-1] * (ZLPTS[0] - Zaxis) + RLPTS[0] * (Zaxis - ZLPTS[NLPTS-1])) /(ZLPTS[0] - ZLPTS[NLPTS-1]);
      else if (RLPTS[NLPTS-1] > Raxis)
	ROUT = (RLPTS[NLPTS-1] * (ZLPTS[0] - Zaxis) + RLPTS[0] * (Zaxis - ZLPTS[NLPTS-1])) /(ZLPTS[0] - ZLPTS[NLPTS-1]);
    }
  RIN  *= R0;
  ROUT *= R0;
  printf ("R0      = %11.4e  B0          = %11.4e\n", R0,  B0);
  printf ("RIN     = %11.4e  ROUT        = %11.4e\n", RIN, ROUT);

  // ........
  // Read Psi
  // ........
  PSIARRAY.resize (NRPTS, NZPTS);  // Psi (R, Z)

  double val; int ival;
  file = OpenFiler ((char*) "Outputs/Stage1/PsiSequential.txt");
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {	
	if (fscanf (file, "%d %d %lf", &ival, &ival, &val) != 3)
	  {
	    printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/PsiSequential.txt\n");
	    exit (1);
	  }
	PSIARRAY (i, j) = val;
      }
  fclose (file);

  // ....................................................
  // Modify PSI such that PsiBoundary = 0 and PsiAxis = 1
  // ....................................................
  double PsiAxis     = InterpolatePsi (Raxis,    Zaxis,    0);
  double PsiBoundary = InterpolatePsi (RBPTS[0], ZBPTS[0], 0);

  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {
	double val = PSIARRAY (i, j) - PsiBoundary;
	PSIARRAY (i, j) = val;
      }

  PsiAxis     = InterpolatePsi (Raxis,    Zaxis,    0);
  PsiBoundary = InterpolatePsi (RBPTS[0], ZBPTS[0], 0);
  Psic        = PsiAxis / (R0*R0*B0);

  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {
	double val = PSIARRAY (i, j) /PsiAxis;
	PSIARRAY (i, j) = val;
      }

  printf ("PsiAxis = %11.4e  PsiBoundary = %11.4e  PsiAxis /(R0*R0*B0) = %11.4e\n", PsiAxis, PsiBoundary, Psic);

  // ......................
  // Output Psi, PsiR, PsiZ
  // ......................
  file = OpenFilew ((char*) "Outputs/Stage2/Psi.txt");
  for (int i = 0; i < NRPTS; i++)
    {
      for (int j = 0; j < NZPTS; j++)
	fprintf (file, "%16.9e ",  PSIARRAY (i, j));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFilew ((char*) "Outputs/Stage2/PsiR.txt");
  for (int i = 0; i < NRPTS; i++)
    {
      for (int j = 0; j < NZPTS; j++)
	fprintf (file, "%16.9e ", GetPsiR (RPTS[i], ZPTS[j]));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFilew ((char*) "Outputs/Stage2/PsiZ.txt");
  for (int i = 0; i < NRPTS; i++)
    {
      for (int j = 0; j < NZPTS; j++)
	fprintf (file, "%16.9e ", GetPsiZ (RPTS[i], ZPTS[j]));
      fprintf (file, "\n");
    }
  fclose (file);

  // .............................
  // Read equilibrium profile data
  // .............................
  PSIN = new double[NRPTS]; // PSI_N array 
  G    = new double[NRPTS]; // g
  Pr   = new double[NRPTS]; // p
  GGp  = new double[NRPTS]; // g dg/dpsi
  Prp  = new double[NRPTS]; // dp/dpsi
  Q    = new double[NRPTS]; // q

  file = OpenFiler ((char*) "Outputs/Stage1/Profiles.txt");
  for (int i = 0; i < NRPTS; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf", &PSIN[i], &G[i], &Pr[i], &GGp[i], &Prp[i], &Q[i]) != 6)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Profiles.txt\n");
	exit (1);
      }
  fclose (file);

  // ....................................................
  // Find closest equilibrium grid-point to magnetic axis
  // ....................................................
  double rmin = 1.e6;
  for (int i = 0; i < NRPTS; i++)
    if (fabs (RPTS[i] - Raxis) < rmin)
      {
	rmin = fabs (RPTS[i] - Raxis);
	ic   = i;
      }
  rmin = 1.e6;
  for (int j = 0; j < NZPTS; j++)
    if (fabs (ZPTS[j] - Zaxis) < rmin)
      {
	rmin = fabs (ZPTS[j] - Zaxis);
	jc   = j;
      }

  // ..........................................................................
  // Find closest equilibrium grid-point to inner magnetic boundary on midplane
  // ..........................................................................
  ia = 0;
  for (int i = 0; i <= ic; i++)
    {
      if (PSIARRAY (i, jc) < 0.)
	ia++;
    }
  Rbound = (RPTS[ia-1]*PSIARRAY (ia, jc) - RPTS[ia]*PSIARRAY (ia-1, jc))
    /(PSIARRAY (ia, jc) - PSIARRAY (ia-1, jc));

  L = ic-ia+2;  // Number of points in Psi (R, Zaxis) array

  // ..........................................................................
  // Find closest equilibrium grid-point to outer magnetic boundary on midplane
  // ..........................................................................
  int ia1 = ic;
  for (int i = ic; i < NRPTS; i++)
    {
      if (PSIARRAY (i, jc) > 0.)
	ia1++;
    }
  Rbound1 = (RPTS[ia1+1]*PSIARRAY (ia1, jc) - RPTS[ia1]*PSIARRAY (ia1+1, jc))
    /(PSIARRAY (ia1, jc) - PSIARRAY (ia1+1, jc));
}

// ######################################
// Function to calculate Stage2 q profile
// ######################################
void Flux::Stage2CalcQ ()
{
  // ............................
  // Setup Psi (R, Zaxis) profile
  // ............................
  s = new double[L]; // Array of s = sqrt [1 - Psi (R, Zaxis)] values

  for (int l = 0; l < L-2; l++)
    s[L-2-l] = sqrt (1. - PSIARRAY (l+ia, jc));
  s[0]   = 0.;
  s[L-1] = 1.;

  // .........................
  // Setup R (Z=Zaxis) profile
  // .........................
  Rs = new double[L]; // Array of R(s) values

  for (int l = 0; l < L-2; l++)
    Rs[L-2-l] = RPTS[l+ia];
  Rs[0]   = Raxis;
  Rs[L-1] = Rbound;

  FILE* file = OpenFilew ((char*) "Outputs/Stage2/rs.txt");
  for (int l = 0; l  < L; l++)
    fprintf (file, "%16.9e %16.9e\n", Rs[l], s[l]);
  fclose (file);

  // ......................
  // Set up Stage2 Psi grid
  // ......................
  P   = new double[NPSI];  // Psi array
  RP  = new double[NPSI];  // R(Psi) on inboard midplane
  RP1 = new double[NPSI];  // R(Psi) on outboard midplane
  Bt  = new double[NPSI];  // B_toroidal(Psi) on inboard midplane
  Bt1 = new double[NPSI];  // B_toroidal(Psi) on outboard midplane
  Bp  = new double[NPSI];  // B_poloidal(Psi) on inboard midplane
  Bp1 = new double[NPSI];  // B_poloidal(Psi) on outboard midplane
  rP  = new double[NPSI];  // r(Psi)
  GP  = new double[NPSI];  // g(Psi)
  QGP = new double[NPSI];  // q(Psi)/g(Psi) 
  QP  = new double[NPSI];  // q(Psi)
  PP  = new double[NPSI];  // P(Psi)
  GPP = new double[NPSI];  // dg/dPsi
  PPP = new double[NPSI];  // dP/dPsi
  S   = new double[NPSI];  // sqrt (1 - Psi)
  QX  = new double[NPSI];  // q(Psi) from gFile

  PsiN = new double[NPSI]; // PsiN array
  QPN  = new double[NPSI]; // dQ/dPsiN array
  A1   = new double[NPSI]; // QP/QPN/fabs(Psic) array
  
  for (int j = 0; j < NPSI; j++)
    {
      double s = double (j) /double (NPSI-1);

      PsiN[j] = 1. - pow (1. - s, PACK);
      P[j]    = 1. - PsiN[j];
      S[j]    = sqrt (1. - P[j]);
    }

  // .......................................
  // Calculate Stage2 g(Psi), P(Psi) profile
  // .......................................
  for (int j = 0; j < NPSI; j++)
    {
      double pval = 1. - P[j];

      GP [j] = Interpolate (NRPTS, PSIN, G,   pval, 0);
      PP [j] = Interpolate (NRPTS, PSIN, Pr,  pval, 0);
      GPP[j] = Interpolate (NRPTS, PSIN, GGp, pval, 0) /GP[j];
      PPP[j] = Interpolate (NRPTS, PSIN, Prp, pval, 0);
      QX [j] = Interpolate (NRPTS, PSIN, Q,   pval, 0);
    }

  // ...............................
  // Calculate Stage2 R(Psi) profile
  // ...............................
  for (int j = 0; j < NPSI; j++)
    RP[j] = Interpolate (L, s, Rs, S[j], 0);
  
  // ......................................
  // Calculate Stage2 q(Psi)/g(Psi) profile
  // ......................................
  printf ("Calculating q(Psi)/g(Psi) profile:\n");
  CalcQGP ();
  
  RP1[0]      = Raxis;
  RP1[NPSI-1] = Rbound1;

  QGP[0] = Q[0] /G[0];
  QP [0] = Q[0];

  // Extrapolate q profile
  QGP[NPSI-1] = InterpolateQ (NPSI-1, PsiN, QGP, 1., 0, 4);
  QP [NPSI-1] = QGP[NPSI-1] * G[NRPTS-1];

  // Calculate and smooth QPN profile (quartic interpolation)
  int NSMOOTH = 8;
  double* QPSmooth = new double[NPSI];

  for (int i = 0; i < NPSI; i++)
    QPSmooth[i] = QP[i];

  for (int i = 0; i < NSMOOTH; i++)
    Smoothing (NPSI, QPSmooth);

   for (int j = 0; j < NPSI; j++)
    QPN [j] = InterpolateQ (NPSI, PsiN, QP, PsiN[j], 1, 4);
   
  for (int i = 0; i < NSMOOTH; i++)
    Smoothing (NPSI, QPN);

  delete[] QPSmooth;
  
  // Calculate A1 profile
  for (int j = 0; j < NPSI; j++)
    A1[j] = QP[j] /QPN[j] /fabs(Psic);

  // ...........................................
  // Calculate B_tor, B_pol profiles on midplane
  // ...........................................
  for (int j = 0; j < NPSI; j++)
    {
      Bt[j]  = GP[j] /RP[j];
      Bt1[j] = GP[j] /RP1[j];

      double PsiR = GetPsiR (RP[j], Zaxis);
      double PsiZ = GetPsiZ (RP[j], Zaxis);

      Bp[j] = fabs(Psic) * sqrt (PsiR*PsiR + PsiZ*PsiZ) /RP[j];

      PsiR = GetPsiR (RP1[j], Zaxis);
      PsiZ = GetPsiZ (RP1[j], Zaxis);

      Bp1[j] = fabs(Psic) * sqrt (PsiR*PsiR + PsiZ*PsiZ) /RP1[j];
    }

  // .............................
  // Calculate Stage2 r(P) profile
  // .............................
  printf ("Calculating r(Psi) profile:\n");
  CalcrP ();
  rP[0] = 0.;
  ra    = rP[NPSI-1];

  // Perform diagnostic integration check (passes!)
  //for (int j = 1; j < NPSI; j++)
  //  printf ("r/a = %11.4e  PsiN = %11.4e  dPsiN/dr = %11.4e  r*g/|psi_c|/q = %11.4e  ratio = %11.4e\n",
  //	    rP[j]/ra, PsiN[j], Interpolate (NPSI, rP, PsiN, rP[j], 1), rP[j]*GP[j]/fabs(Psic)/QP[j], rP[j]*GP[j]/fabs(Psic)/QP[j] /Interpolate (NPSI, rP, PsiN, rP[j], 1));

  // .........................
  // Find q95, r95, qlim, rlim
  // .........................
  double g95;
  double s95  = sqrt (0.95);
  double slim = sqrt (PSILIM);
  double sa   = 1.;
  q95    = Interpolate (NPSI, S, QP, s95,  0);
  r95    = Interpolate (NPSI, S, rP, s95,  0);
  g95    = Interpolate (NPSI, S, GP, s95,  0);
  qlim   = Interpolate (NPSI, S, QP, slim, 0);
  rlim   = Interpolate (NPSI, S, rP, slim, 0);
  qa     = Interpolate (NPSI, S, QP, sa,   0);
  printf ("q95 = %11.4e  r95/ra = %11.4e  qlim = %11.4e  rlim/ra = %11.4e  qa = %11.4e  ra = %11.4e  a = %11.4e (m)\n",
	  q95, r95 /ra, qlim, rlim/ra, qa, ra, ra*R0);

  // ..................
  // Output q95 and r95
  // ..................
  file = OpenFilew ((char*) "Outputs/Stage2/q95.txt");
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", q95, r95 /ra, QP[0], QP[NPSI-1], qlim, rlim /ra);
  fclose (file);

  // ......................
  // Output Stage2 profiles
  // ......................
  file = OpenFilew ((char*) "Outputs/Stage2/qr.txt");
  for (int j = 0; j < NPSI; j++)
    {
      if (1. - P[j] < PSILIM)
	fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", 
		 rP[j] /ra, QP[j], QGP[j], P[j], GP[j], PP[j], GPP[j], PPP[j], QX[j], RP[j], RP1[j], Bt[j], Bt1[j], Bp[j], Bp1[j], QPN[j], A1[j]);
    }
  fclose (file);
}

// ###################################
// Function to find rational surfaces.
// Assumes monotonic q-profile.
// ###################################
void Flux::Stage2FindRational ()
{
  // .....................................
  // Determine number of rational surfaces
  // .....................................
  int mmin = int (double (NTOR) * QP[0]) + 1;
  int mmax = int (double (NTOR) * QP[NPSI-1]);
  if (mmin < MMIN)
    mmin = MMIN;
  if (mmax > MMAX)
    mmax = MMAX;
  nres = mmax - mmin + 1;

  if (nres <= 0)
    {
      printf ("FLUX::Stage2FindRational: Error no rational surfaces in plasma\n");
      exit (1);
    }

  // .....................................
  // Calculate rational surface quantities
  // .....................................
  mres    = new int   [nres];
  qres    = new double[nres];
  rres    = new double[nres];
  sres    = new double[nres];
  gres    = new double[nres];
  Rres    = new double[nres];
  Rres1   = new double[nres];
  gmres   = new double[nres];
  fcres   = new double[nres];
  PsiNres = new double[nres];
  A1res   = new double[nres];
  for (int i = 0; i < nres; i++)
    {
      mres[i] = mmin + i;

      qres[i] = double (mres[i]) /double (NTOR);
      rres[i] = Interpolate (NPSI, QP, rP, qres[i], 0);
      sres[i] = rres[i] /Interpolate (NPSI, QP, rP, qres[i], 1) /qres[i];

      // Check that PsiN < PSILIM
      PsiNres[i] = 1. - Interpolate (NPSI, rP, P, rres[i], 0);
      if (PsiNres[i] > PSILIM)
	{
	  nres = i;
	  break; 
	}

      // Correct rres values
      for (int ii = 0; ii < 4; ii++)
	{
	  double qqq = Interpolate (NPSI, rP, QP, rres[i], 0);
	  double qqp = Interpolate (NPSI, rP, QP, rres[i], 1);
	  double qpp = Interpolate (NPSI, rP, QP, rres[i], 2);
	  rres[i] += (- qqp + sqrt (qqp*qqp - 2.*qpp * (qqq - qres[i]))) /qpp;
	}
      
      gmres  [i] = Interpolate (NPSI, rP, QP, rres[i], 0) - qres[i];
      sres   [i] = rres[i] /Interpolate (NPSI, QP, rP, qres[i], 1) /qres[i];
      gres   [i] = Interpolate (NPSI, rP, GP,   rres[i], 0);
      Rres   [i] = Interpolate (NPSI, rP, RP,   rres[i], 0);
      Rres1  [i] = Interpolate (NPSI, rP, RP1,  rres[i], 0);
      PsiNres[i] = Interpolate (NPSI, rP, PsiN, rres[i], 0);
      A1res  [i] = Interpolate (NPSI, rP, A1,   rres[i], 0);
    }

  printf ("Rational surface data:\n");
  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d  PsiNs = %11.4e  rs/ra = %11.4e  ss = %11.4e  residual = %11.4e  R = %11.4e  A1 = %11.4e  residual = %11.4e\n",
	    mres[i], PsiNres[i], rres[i] /ra, sres[i], gmres[i], Rres1[i], A1res[i], 1. - (gres[i]/sres[i]/qres[i]) *rres[i]*rres[i] /Psic/Psic /A1res[i]);
   
  // .....................................
  // Confirm q values at rational surfaces
  // .....................................
  printf ("Confirm rational surface q-values:\n");
  CheckQP ();

  // ................................
  // Find appropriate PSILIM for GPEC
  // ................................
  double dmlim = 0.2;
  double qGPEC = qres[nres-1] + dmlim /double (NTOR);
  double rGPEC = Interpolate (NPSI, QP, rP, qGPEC, 0);
  for (int ii = 0; ii < 4; ii++)
	{
	  double qqq = Interpolate (NPSI, rP, QP, rGPEC, 0);
	  double qqp = Interpolate (NPSI, rP, QP, rGPEC, 1);
	  double qpp = Interpolate (NPSI, rP, QP, rGPEC, 2);
	  rGPEC += (- qqp + sqrt (qqp*qqp - 2.*qpp * (qqq - qGPEC))) /qpp;
	}
  double PSIGPEC = Interpolate (NPSI, rP, PsiN, rGPEC, 0);
  double qqGPEC  = Interpolate (NPSI, rP, QP,   rGPEC, 0);
  printf ("Calculating PSILIM for GPEC:\n");
  printf ("qGPEC = %11.4e  PSIGPEC = %11.4e  res = %11.4e\n", qGPEC, PSIGPEC, fabs (qqGPEC - qGPEC));
 }

// ##################################################
// Calculate straight angle data at rational surfaces
// ##################################################
void Flux::Stage2CalcStraightAngle ()
{
  // ..................
  // Set up theta array
  // ..................
  th = new double[NTHETA]; 
  for (int k = 0; k < NTHETA; k++)
    {
      double t = double (k) /double (NTHETA-1);
      th[k]    = 2.* M_PI * t;
    }
  Rst = gsl_matrix_alloc (nres, NTHETA);
  Zst = gsl_matrix_alloc (nres, NTHETA);

  // .............................
  // Calculate straight angle data
  // .............................
  printf ("Calculating straight angle data at rational surface:\n");
  CalcStraightAngle ();

  // ...........
  // Output nres
  // ...........
  FILE* file = OpenFilew ((char*) "Outputs/Stage2/nres.txt");
  fprintf (file, "%d\n", nres);
  fclose (file);

  // ..........
  // Output Rst
  // ..........
  file = OpenFilew ((char*) "Outputs/Stage2/Rst.txt");
  for (int k = 0; k < NTHETA; k++)
    {
      fprintf (file, "%16.9e ", th[k] /M_PI);
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_matrix_get (Rst, j, k));
      fprintf (file, "\n");
    }
  fclose (file);

  // ..........
  // Output Zst
  // ..........
  file = OpenFilew ((char*) "Outputs/Stage2/Zst.txt");
  for (int k = 0; k < NTHETA; k++)
    {
      fprintf (file, "%16.9e ", th[k] /M_PI);
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_matrix_get (Zst, j, k));
      fprintf (file, "\n");
    }
  fclose (file);
}

// ######################################################
// Calculate neoclassical angle data at rational surfaces
// ######################################################
void Flux::Stage2CalcNeoclassicalAngle ()
{
  // ...........................................
  // Calculate gamma values at rational surfaces
  // ...........................................
  printf ("Calculating gamma values at rational surfaces:\n");
  CalcGamma ();

  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d  rs/ra = %11.4e  g = %11.4e  gamma = %11.4e  gamma*q/g = %11.4e\n",
	    mres[i], rres[i] /ra, gres[i], gmres[i], gmres[i]*qres[i] /gres[i]);

  // ..................
  // Set up Theta array
  // ..................
  Th = new double[NTHETA];
  for (int k = 0; k < NTHETA; k++)
    {
      double t = double (k) /double (NTHETA-1);
      Th[k]    = 2.* M_PI * t;
    }
  Rnc = gsl_matrix_alloc (nres, NTHETA);
  Znc = gsl_matrix_alloc (nres, NTHETA);

  // .................................
  // Calculate neoclassical angle data
  // .................................
  printf ("Calculating neoclassical angle data at rational surfaces:\n");
  CalcNeoclassicalAngle ();

  // ..........
  // Output Rnc
  // ..........
  FILE* file = OpenFilew ((char*) "Outputs/Stage2/Rnc.txt");
  for (int k = 0; k < NTHETA; k++)
    {
      fprintf (file, "%16.9e ", Th[k] /M_PI);
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_matrix_get (Rnc, j, k));
      fprintf (file, "\n");
    }
  fclose (file);

  // ..........
  // Output Znc
  // ..........
  file = OpenFilew ((char*) "Outputs/Stage2/Znc.txt");
  for (int k = 0; k < NTHETA; k++)
    {
      fprintf (file, "%16.9e ", Th[k] /M_PI);
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_matrix_get (Znc, j, k));
      fprintf (file, "\n");
    }
  fclose (file);

  // ..................................
  // Calculate |B| on rational surfaces
  // ..................................
  Bnc = gsl_matrix_alloc (nres, NTHETA);
  for (int j = 0; j < nres; j++)
    for (int k = 0; k < NTHETA; k++)
      {
	double gval = gres[j];
	double Rval = gsl_matrix_get (Rnc, j, k);
	double Zval = gsl_matrix_get (Znc, j, k);
	double PsiR = GetPsiR (Rval, Zval);
	double PsiZ = GetPsiZ (Rval, Zval);
	double Grad = sqrt (PsiR*PsiR + PsiZ*PsiZ);
	double Bval = sqrt (gval*gval + Psic*Psic * Grad*Grad) /Rval;

	gsl_matrix_set (Bnc, j, k, Bval);
      }

  // ..........
  // Output Bnc
  // ..........
  file = OpenFilew ((char*) "Outputs/Stage2/Bnc.txt");
  for (int k = 0; k < NTHETA; k++)
    {
      fprintf (file, "%16.9e ", Th[k] /M_PI);
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_matrix_get (Bnc, j, k));
      fprintf (file, "\n");
    }
  fclose (file);

  // .........................................
  // Calculate d|B|dTheta on rational surfaces
  // .........................................
  Cnc = gsl_matrix_alloc (nres, NTHETA);
  double* Y = new double[NTHETA];
  for (int j = 0; j < nres; j++)
    {
      for (int k = 0; k < NTHETA; k++)
 	Y[k] = gsl_matrix_get (Bnc, j, k);

      for (int k = 0; k < NTHETA; k++)
	{
	  double val = InterpolatePeriodic (NTHETA, Th, Y, Th[k], 1);
	  gsl_matrix_set (Cnc, j, k, val);
	}
    }
  delete[] Y;

  // ..........
  // Output Cnc
  // ..........
  file = OpenFilew ((char*) "Outputs/Stage2/Cnc.txt");
  for (int k = 0; k < NTHETA; k++)
    {
      fprintf (file, "%16.9e ", Th[k]/M_PI);
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_matrix_get (Cnc, j, k));
      fprintf (file, "\n");
    }
  fclose (file);
}

// ######################################################
// Calculate neoclassical parameters at rational surfaces
// ######################################################
void Flux::Stage2CalcNeoclassicalPara ()
{
  // ................
  // Allocate memory
  // ...............
  I1     = new double[nres];
  I2     = new double[nres];
  I3     = new double[nres];
  I4     = gsl_matrix_alloc (nres, NNC);
  I5     = gsl_matrix_alloc (nres, NNC);
  I6     = gsl_matrix_alloc (nres, NTHETA);
  Ktres  = new double[nres];
  Kares  = new double[nres];
  Khres  = new double[nres];
  ajj    = new double[nres];
  dPsidr = new double[nres];

  // ......................................................
  // Calculate neoclassical parameters at rational surfaces
  // ......................................................
  printf ("Calculating neoclassical parameters at rational surfaces:\n");
  double sum;
  for (int i = 0; i < nres; i++)
    {
      // Calculate I1
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) /gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I1[i] = sum;

      // Calculate I2
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) * gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I2[i] = sum;

      // Calculate I3
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) * gsl_matrix_get (Cnc, i, j) * gsl_matrix_get (Cnc, i, j)
	  /gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I3[i] = sum;

      // Calculate I4
      for (int k = 0; k < NNC; k++)
	{
	  double kk = double (k+1);

	  sum = 0.;
	  for (int j = 0; j < NTHETA; j++)
	    sum += Weight1D (j) * cos (kk * Th[j]) /gsl_matrix_get (Bnc, i, j);
	  sum *= sqrt (2.*kk) /2./M_PI;

	  gsl_matrix_set (I4, i, k, sum);
	}

      // Calculate I5
      for (int k = 0; k < NNC; k++)
	{
	  double kk = double (k+1);

	  sum = 0.;
	  for (int j = 0; j < NTHETA; j++)
	    sum += Weight1D (j) * cos (kk * Th[j]) /gsl_matrix_get (Bnc, i, j)
	      /gsl_matrix_get (Bnc, i, j) /2.;
	  sum *= sqrt (2.*kk) /2./M_PI;

	  gsl_matrix_set (I5, i, k, sum);
	}

      // Calcuate Kt, Ka, and Kh
      sum = 0.;
      for (int k = 0; k < NNC; k++)
	sum += gsl_matrix_get (I4, i, k) * gsl_matrix_get (I5, i, k);
      Ktres[i] = I1[i]*I1[i] * I3[i] /I2[i]/I2[i] /sum;
      Kares[i] = (8./3./M_PI) * (I2[i] /I3[i]) * Ktres[i]*Ktres[i];
      Khres[i] = gres[i]*gres[i] * I1[i] /Rres1[i]/Rres1[i] /I2[i];

      // Calculate Bmax
      double Bmax = -1.;
      for (int k = 0; k < NTHETA; k++)
	if (gsl_matrix_get (Bnc, i, k) > Bmax)
	  Bmax = gsl_matrix_get (Bnc, i, k);
      
      // Calculate I6
      for (int k = 0; k < NTHETA; k++)
	{
	  double lambda = double (k) * hh;

	  sum = 0.;
	  for (int kk = 0; kk < NTHETA; kk++)
	    sum += Weight1D (kk) * sqrt (1. - lambda * gsl_matrix_get (Bnc, i, kk) /Bmax)
	      /gsl_matrix_get (Bnc, i, kk);
	  sum /= 2.*M_PI;
	  
	  gsl_matrix_set (I6, i, k, sum);
	}

      // Calculate fraction of circulating particles
      sum = 0.;
      for (int k = 0; k < NTHETA; k++)
	{
	  double lambda = double (k) * hh;
	  sum += weight1D (k) * lambda /gsl_matrix_get (I6, i, k);
	}
      sum *= 0.75 * I2[i] /Bmax/Bmax;

      fcres[i] = sum;

      // Calculate ajj values
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	{
	  double Rval = gsl_matrix_get (Rst, i, j);
	  double Zval = gsl_matrix_get (Zst, i, j);
 	  double PsiR = GetPsiR (Rval, Zval);
	  double PsiZ = GetPsiZ (Rval, Zval);
	  sum += Weight1D (j) /(PsiR*PsiR + PsiZ*PsiZ);
	}
      double fac = rres[i] * gres[i] /qres[i] /Psic;
      sum *= fac*fac /(2.*M_PI);
     
      ajj[i] = sum;

      // Calculate dPsidr values
      dPsidr[i] = rres[i] * gres[i] /qres[i] /fabs(Psic);
    }

  // ..............................
  // Output neoclassical parameters
  // ..............................
  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d I1 = %11.4e I2 = %11.4e I3 = %11.4e I4 = (%11.4e; %11.4e) I5 = (%11.4e; %11.4e) Kt = %11.4e Ka = %11.4e Kh = %11.4e fc = %11.4e ajj = %11.4e dPsiNdr = %11.4e\n",
	    mres[i], I1[i], I2[i], I3[i],
	    gsl_matrix_get (I4, i, 0), gsl_matrix_get (I4, i, NNC-1),
	    gsl_matrix_get (I5, i, 0), gsl_matrix_get (I5, i, NNC-1),
	    Ktres[i], Kares[i], Khres[i], fcres[i], ajj[i], dPsidr[i]);
}

// ############################
// Calculate stability matrices
// ############################
void Flux::Stage2CalcMatrices ()
{
  printf ("Calculating stability matrices:\n");
  
  // ..................
  // Calculate F matrix
  // ..................
  FF = gsl_matrix_complex_alloc (nres, nres);
  EE = gsl_matrix_complex_alloc (nres, nres);
 
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      {
 	double sumc = 0.;
	for (int k = 0; k < NTHETA; k++)
	  for (int kk = 0; kk < NTHETA; kk++)
	    sumc += Weight2D (k, kk) * GreenPlasmaCos (i, k, j, kk);
	sumc /= 2.*M_PI * 2.*M_PI;

	double sums = 0.;
	for (int k = 0; k < NTHETA; k++)
	  for (int kk = 0; kk < NTHETA; kk++)
	    sums += Weight2D (k, kk) * GreenPlasmaSin (i, k, j, kk);
	sums /= 2.*M_PI * 2.*M_PI;

	gsl_matrix_complex_set (FF, i, j, gsl_complex_rect (sumc, sums));
      }

  printf ("F-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%9.2e,%9.2e) ", GSL_REAL (gsl_matrix_complex_get (FF, i, j)), GSL_IMAG (gsl_matrix_complex_get (FF, i, j)));
      printf ("\n");
    }

  FILE* file = OpenFilew ("Outputs/Stage2/F_Matrix.txt");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_complex_abs (gsl_matrix_complex_get (FF, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);

  // ------------------
  // Calculate E-matrix
  // ------------------
  gsl_permutation*    px  = gsl_permutation_alloc (nres);
  int                 sss = 0;
  gsl_matrix_complex* FFF = gsl_matrix_complex_alloc (nres, nres);

  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      gsl_matrix_complex_set (FFF, i, j, gsl_matrix_complex_get (FF, i, j));

  gsl_linalg_complex_LU_decomp (FFF, px, &sss);
  gsl_linalg_complex_LU_invert (FFF, px, EE);

  gsl_permutation_free    (px);
  gsl_matrix_complex_free (FFF);

  printf ("E-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%9.2e,%9.2e) ", GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
      printf ("\n");
    }

  file = OpenFilew ("Outputs/Stage2/E_Matrix.txt");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	fprintf (file, "%16.9e ", gsl_complex_abs (gsl_matrix_complex_get (EE, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);

  // -------------------
  // Calculate E-vectors
  // -------------------
  EI = gsl_vector_complex_alloc (nres);
  EO = gsl_vector_complex_alloc (nres);
 
  for (int i = 0; i < nres; i++)
    {
      double sumc = 0.;
      for (int j = 0; j < NTHETA; j++)
	sumc += Weight1D (j) * GreenInboardCos (i, j);
      sumc *= (R0 /RIN) /2./M_PI;

      double sums = 0.;
      for (int j = 0; j < NTHETA; j++)
	sums += Weight1D (j) * GreenInboardSin (i, j);
      sums *= (R0 /RIN) /2./M_PI;

      gsl_vector_complex_set (EI, i, gsl_complex_rect (sumc, sums));
   
      sumc = 0.;
      for (int j = 0; j < NTHETA; j++)
	sumc += Weight1D (j) * GreenOutboardCos (i, j);
      sumc *= (R0 /ROUT) /2./M_PI;

      sums = 0.;
      for (int j = 0; j < NTHETA; j++)
	sums += Weight1D (j) * GreenOutboardSin (i, j);
      sums *= (R0 /ROUT) /2./M_PI;

      gsl_vector_complex_set (EO, i, gsl_complex_rect (sumc, sums));
    }

  printf ("E-vectors:\n");
  for (int i = 0; i < nres; i++)
    printf ("mpol = %4d  rs/ra = %11.4e  EI = (%11.4e, %11.4e)  EO = (%11.4e, %11.4e)\n", mres[i], rres[i]/ra,
	    GSL_REAL (gsl_vector_complex_get (EI, i)), GSL_IMAG (gsl_vector_complex_get (EI, i)),
	    GSL_REAL (gsl_vector_complex_get (EO, i)), GSL_IMAG (gsl_vector_complex_get (EO, i)));
}
