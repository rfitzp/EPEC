// Neoclassical.cpp

// PROGRAM ORGANIZATION:
// 
//        Neoclassical:: Neoclassical         ()
// void   Neoclassical:: Solve                ()
// void   Neoclassical:: Read_Parameters      ()
// void   Neoclassical:: Read_Profiles        ()
// void   Neoclassical:: Get_Derived          ()
// void   Neoclassical:: Get_Viscosities      ()
// void   Neoclassical:: Get_Parameters       ()
// void   Neoclassical:: Get_Frequencies      ()
// void   Neoclassical:: Get_LayerWidths      ()
// void   Neoclassical:: WriteStage3Netcdfc   ()
// void   Neoclassical:: Get_Normalized       ()
// double Neoclassical:: psi_fun              (double x)
// double Neoclassical:: psi_fun_p            (double x)
// void   Neoclassical:: RK4Adaptive          (double& x, Array<double,1>& y, double& h, 
//			                       double& t_err, double acc, double S, int& rept,
//			                       int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
// void   Neoclassical:: RK4Fixed             (double& x, Array<double,1>& y, double h)
// void   Neoclassical:: Matrix_Mult          (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i)
// void   Neoclassical:: Matrix_Add           (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i)
// void   Neoclassical:: Matrix_Sub           (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i)
// FILE*  Neoclassical:: OpenFilew            (char* filename)
// FILE*  Neoclassical:: OpenFilew            (char* filename)
// FILE*  Neoclassical:: OpenFilea            (char* filename)

#include "Neoclassical.h"

// ###########
// Constructor
// ###########
Neoclassical::Neoclassical ()
{ 
  // ----------------------
  // Set physical constants
  // ----------------------
  e         = GSL_CONST_MKSA_ELECTRON_CHARGE;
  epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
  mu_0      = GSL_CONST_MKSA_VACUUM_PERMEABILITY;
  m_p       = GSL_CONST_MKSA_MASS_PROTON;
  m_e       = GSL_CONST_MKSA_MASS_ELECTRON;

  // -----------------------------------------------------
  // Set default values of adaptive integration parameters
  // -----------------------------------------------------
  xmin    = 1.e-7;
  xmax    = 100.;
  h0      = 1.e-4;
  acc     = 1.e-11;
  hmin    = 1.e-10;
  hmax    = 1.e-2;
  maxrept = 50;
}

// ##############
// Solve problem
// ##############
void Neoclassical::Solve ()
{
  // Read discharge parameters
  Read_Parameters ();
  fflush (stdout);
  
  // Read equilibrium data
  Read_Equilibrium ();
  fflush (stdout);

  // Read profile data
  Read_Profiles ();
  fflush (stdout);

  // Calculate derived quantities at rational surface
  Get_Derived ();
  fflush (stdout);

  // Calculate neoclassical viscosities at rational surfaces
  Get_Viscosities ();
  fflush (stdout);

  // Calculate neoclassical parameters at rational surfaces
  Get_Parameters ();
  fflush (stdout);

  // Calculate neoclassical frequencies at rational surfaces
  Get_Frequencies ();
  fflush (stdout);

  // Calculate linear layer widths at rational surfaces
  Get_LayerWidths ();
  fflush (stdout);

  // Output NETCDF file
  WriteStage3Netcdfc ();

  // Calculate normalized quantites at rational surface and output nFile
  Get_Normalized ();
  fflush (stdout);
}

// ####################################
// Read Neoclassical control parameters
// ####################################
void Neoclassical::Read_Parameters ()
{
  // Set default values of control parameters
  EXB      = 0;

  INTP     = 0;
  INTF     = 0;
  INTC     = 0;

  IMPURITY = 1;

  NEUTRAL  = 1;
  NTYPE    = 0;
  NN       = 0.;
  LN       = 1.;
  SVN      = 0.;
  YN       = 1.;
  EN       = 0.;

  CATS     = 0;

  TIME     = 0.;

  COULOMB  = 17.;
  NSMOOTH  = 100;

  TAUMIN   = - 0.8;
  DMAX     = 10.;
  DMIN     = 0.01;

  WPSIMAX  = 0.2;
 
  // Read namelist
  NameListRead (&IMPURITY, &NEUTRAL, &EXB, &INTP, &INTF, &INTC, &NTYPE, &NN, &LN, &SVN, &YN, &EN, &TIME,
		&COULOMB, &NSMOOTH, &CATS, &TAUMIN, &DMIN, &DMAX, &WPSIMAX);

  // Sanity check
  if (YN < 0.)
    {
      printf ("NEOCLASSICAL::Read_Parameters: Error YN must be positive\n");
      exit (1);
    }
  if (NN < 0.)
     {
       printf ("NEOCLASSICAL::Read_Parameters: Error NN must be positive\n");
       exit (1);
     }
  if (LN < 0.)
     {
       printf ("NEOCLASSICAL::Read_Parameters: Error LN must be positive\n");
       exit (1);
     }
  if (SVN < 0.)
     {
       printf ("NEOCLASSICAL::Read_Parameters: Error SVN must be positive\n");
      exit (1);
     }
  if (EN < 0.)
     {
       printf ("NEOCLASSICAL::Read_Parameters: Error EN must be positive\n");
      exit (1);
     }
  if (NTYPE < 0 || NTYPE > 1)
     {
       printf ("NEOCLASSICAL::Read_Parameters: Error invalid NTYPE value\n");
       exit (1);
     }
  if (EXB < 0 || EXB > 2)
    {
      printf ("NEOCLASSICAL::Read_Parameters: Error invalid EXB value\n");
      exit (1);
    }
  if (TAUMIN <= -1.)
    {
      printf ("NEOCLASSICAL::Read_Parameters: Error invalid TAUMIN value\n");
      exit (1);
    }
  if (DMIN < 0. || DMAX < DMIN)
     {
      printf ("NEOCLASSICAL::Read_Parameters: Error invalid DMIN/DMAX value\n");
      exit (1);
    }
  if (WPSIMAX <= 0.)
    {
      printf ("NEOCLASSICAL::Read_Parameters: Error invalid WPSIMAX value\n");
      exit (1);
    }

  // Output input parameters
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Input parameters (from Inputs/Neoclassical.nml):\n");
  printf ("IMPURITY = %2d NEUTRAL = %2d EXB = %2d INTP = %2d INTF = %2d INTC = %2d NTYPE = %2d NN = %11.4e LN = %11.4e SVN = %11.4e YN = %11.4e EN = %11.4e TIME = %11.4e NSMOOTH = %3d CATS = %2d TAUMIN = %11.4e DMIN = %11.4e DMAX = %11.4e WPSIMAX = %11.4e\n",
	  IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, NSMOOTH, CATS, TAUMIN, DMIN, DMAX, WPSIMAX);
  
  FILE* namelist = OpenFilew ((char*) "Outputs/InputParameters.txt");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "Input parameters (from Inputs/Neoclassical.nml):\n");
  fprintf (namelist, "IMPURITY = %2d NEUTRAL = %2d EXB = %2d INTP = %2d INTF = %2d INTC = %2d NTYPE = %2d NN = %11.4e LN = %11.4e SVN = %11.4e YN = %11.4e EN = %11.4e TIME = %11.4e NSMOOTH = %3d CATS = %2d TAUMIN = %11.4e DMIN = %11.4e DMAX = %11.4e WPSIMAX = %11.4e\n",
	   IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, NSMOOTH, CATS, TAUMIN, DMIN, DMAX, WPSIMAX);
  fclose (namelist);
}

// ###################################################
// Read plasma equilibrium data output by program FLUX
// ###################################################
void Neoclassical::Read_Equilibrium ()
{
  // .................
  // Interpolate fFile
  // .................
  if (INTF != 0)
    {
      // Save pwd
      char pwd[MAXFILENAMELENGTH];
      getcwd (pwd, MAXFILENAMELENGTH);

      // Remove fFile
      CallSystem ("rm -rf Inputs/fFile");
      
      // Get fFiles directory
      char fFileDir[MAXFILENAMELENGTH];
      CallSystem ("greadlink -f Inputs/fFiles > fFileDir");
      FILE* ffd = OpenFiler ("fFileDir");
      fscanf (ffd, "%s", fFileDir);
      fclose (ffd);
      CallSystem ("rm fFileDir");
        
      // Read fFile data
      char           Basename[MAXFILENAMELENGTH];
      char           Filename[MAXFILENAMELENGTH];
      char           filename[MAXFILENAMELENGTH];
      vector<string> fFileName;
      double         filetime;
      vector<double> fFileTime;
      int            fFileNumber;
      
      printf ("Reading fFile data:\n");

      chdir (fFileDir);
      getcwd (Basename, MAXFILENAMELENGTH);
      strcat (Basename, "/");

      FILE* file = OpenFiler ((char*) "Index");
      
      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  strcpy (Filename, Basename);
	  strcat (Filename, filename);
	  
	  fFileName.push_back (Filename);
	  fFileTime.push_back (filetime);
	}
      fFileNumber = fFileName.size ();

      fclose (file);
      chdir (pwd);
      
      // Interpolate fFiles
      fFileInterp (fFileName, fFileTime, fFileNumber, TIME);
    }

  // Read fFile
  printf  ("Reading fFile:\n");
  
  // Read parameters
  FILE* file = OpenFiler ((char*) "Inputs/fFile");
  double in;
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &R_0, &B_0, &a, &in, &in, &in, &in, &in, &in, &NPSI, &ntor, &nres, &in, &in, &in, &PSIRAT) != 16)
    {
      printf ("NEOCLASSICAL::Read_Equilibrium: Error reading fFile (1)\n");
      exit (1);  
    }
  printf ("ntor = %3d  nres = %3d\n", ntor, nres);
  
  // Read equilibrium profiles
  psi.resize    (NPSI);
  rr.resize     (NPSI);
  dpsidr.resize (NPSI);

  for (int j = 0; j < NPSI; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &psi (j), &rr (j), &dpsidr (j)) != 3)
	{
	  printf ("NEOCLASSICAL::Read_Equilibrium: Error reading fFile (2)\n");
	  exit (1);
	}
    }

  // Read rational surface data
  mk.resize     (nres);
  rk.resize     (nres);
  qk.resize     (nres);
  sk.resize     (nres);
  gk.resize     (nres);
  gmk.resize    (nres);
  Ktk.resize    (nres);
  Kastk.resize  (nres);
  Kthek.resize  (nres);
  fck.resize    (nres);
  akk.resize    (nres);
  PsiNk.resize  (nres);
  dPsidr.resize (nres);
  A2.resize     (nres);
  q_hat.resize  (nres);
  C1.resize     (nres);
  C2.resize     (nres);
  DR.resize     (nres);
  Poem1.resize  (nres);
  Poem2.resize  (nres);
  Poem3.resize  (nres);
  Deltaw.resize (nres);
  Sigmaw.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mk (j), &rk (j), &sk (j), &gk (j), &gmk (j), &Ktk (j), &Kastk (j), &fck (j), &akk (j), &PsiNk (j),
		  &dPsidr (j), &Kthek (j), &in, &A2 (j), &q_hat (j), &C1 (j), &C2 (j), &DR (j), &Poem1 (j),
		  &Poem2 (j), &Poem3 (j), &Deltaw (j), &Sigmaw (j)) != 23)
	{
	  printf ("NEOCLASSICAL:Read_Equilibrium: Error reading fFile (3)\n");
	  exit (1);
	}
      qk (j) = double (mk (j)) /double (ntor);
      printf ("m = %3d  r = %9.2e  s = %9.2e  g = %9.2e  gm = %9.2e  Kt = %9.2e  Kast = %9.2e  Kthe = %9.2e  fc = %9.2e  akk = %9.2e  PsiN = %9.2e  q_hat = %9.2e  DR = %9.2e  Poem1 = %9.2e  Poem2 = %9.2e  Poem3 = %9.2e  Deltaw = %9.2e  Sigmaw = %9.2e\n",
	      mk (j), rk (j), sk (j), gk (j), gmk (j), Ktk (j), Kastk (j), Kthek (j), fck (j), akk (j), PsiNk (j), q_hat (j), DR (j),
	      Poem1 (j), Poem2 (j), Poem3 (j), Deltaw (j), Sigmaw (j));
    }

  int ini;
  for (int j = 0; j < nres; j++)
    for (int k = 0; k < nres; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &in, &in) != 4)
	 {
	   printf ("NEOCLASSICAL::Error reading fFile (4)\n");
	   exit (1);
	 }
      }
  
  EEh.resize (nres, nres);
   
  for (int j = 0; j < nres; j++)
    for (int k = 0; k < nres; k++)
      {
	double Ereal, Eimag;
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal, &Eimag) != 4)
	 {
	   printf ("NEOCLASSCICAL::Error reading fFile (5)\n");
	   exit (1);
	 }

	EEh (j, k) = gsl_complex_abs (gsl_complex_rect (Ereal, Eimag));
      }

  fclose (file);
}

// ############################
// Read profile data from pFile
// ############################
void Neoclassical::Read_Profiles ()
{
  // Interpolate pFiles
  if (INTP != 0)
    {
      // Save pwd
      char pwd[MAXFILENAMELENGTH];
      getcwd (pwd, MAXFILENAMELENGTH);

      // Remove pFile
      CallSystem ("rm -rf Inputs/pFile");

      // Get pFiles directory
      char pFileDir[MAXFILENAMELENGTH];
      CallSystem ("greadlink -f Inputs/pFiles > pFileDir");
      FILE* pfd = OpenFiler ("pFileDir");
      fscanf (pfd, "%s", pFileDir);
      fclose (pfd);
      CallSystem ("rm pFileDir");
  
      // Read pFile data
      char           Basename[MAXFILENAMELENGTH];
      char           Filename[MAXFILENAMELENGTH];
      char           filename[MAXFILENAMELENGTH];
      vector<string> pFileName;
      double         filetime;
      vector<double> pFileTime;
      int            pFileNumber;
  
      printf ("Reading pFile data:\n");

      chdir (pFileDir);
      getcwd (Basename, MAXFILENAMELENGTH);
      strcat (Basename, "/");
   
      FILE* file = OpenFiler ((char*) "Index");
   
      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  strcpy (Filename, Basename);
	  strcat (Filename, filename);
	  
	  pFileName.push_back (Filename);
	  pFileTime.push_back (filetime);
	}
      pFileNumber = pFileTime.size ();

      fclose (file);
      chdir (pwd);

      // Interpolate pFiles
      pFileInterp (pFileName, pFileTime, pFileNumber, TIME);
    }
  
  // Read profile data from pFile
  printf ("Reading profile data:\n");
  pFileRead ();

  // Determine ion data (impurities = 0; majority = 1)
  NI  = int (NZA.GetX (1) + 1.e-6);
  ZI  = NZA.GetY      (1);
  AI  = NZA.GetdYdX   (1);
  NII = int (NZA.GetX (0) + 1.e-6);
  ZII = NZA.GetY      (0);
  AII = NZA.GetdYdX   (0);

  printf ("NI  = %3d  ZI  = %9.2e  AI  = %9.2e\n", NI,  ZI,  AI);
  printf ("NII = %3d  ZII = %9.2e  AII = %9.2e\n", NII, ZII, AII);

  // Interpolate profiles onto existing psi grid
  double fac0 = dpsidr (0)      /a;
  double fac1 = dpsidr (NPSI-1) /a;
  
  n_e.resize    (NPSI); dn_edr.resize (NPSI);
  T_e.resize    (NPSI); dT_edr.resize (NPSI);
  n_i.resize    (NPSI); dn_idr.resize (NPSI);
  T_i.resize    (NPSI); dT_idr.resize (NPSI);
  n_b.resize    (NPSI); w_E.resize    (NPSI);
  w_t.resize    (NPSI); n_I.resize    (NPSI);
  dn_Idr.resize (NPSI); T_I.resize    (NPSI);
  dT_Idr.resize (NPSI); n_n.resize    (NPSI);
  Quasi.resize  (NPSI); Z_eff.resize  (NPSI);
  alpha.resize  (NPSI); w_p.resize    (NPSI);

  dn_edP1.resize (NPSI); dT_edP1.resize (NPSI); dn_idP1.resize (NPSI); dT_idP1.resize (NPSI);
  dn_edP2.resize (NPSI); dT_edP2.resize (NPSI); dn_idP2.resize (NPSI); dT_idP2.resize (NPSI);
  dn_edP3.resize (NPSI); dT_edP3.resize (NPSI); dn_idP3.resize (NPSI); dT_idP3.resize (NPSI);

  n_e    (0)  = ne.GetY    (0);
  dn_edr (0)  = ne.GetdYdX (0) * fac0;
  T_e    (0)  = Te.GetY    (0);
  dT_edr (0)  = Te.GetdYdX (0) * fac0;
  n_i    (0)  = ni.GetY    (0);
  dn_idr (0)  = ni.GetdYdX (0) * fac0;
  T_i    (0)  = Ti.GetY    (0);
  dT_idr (0)  = Ti.GetdYdX (0) * fac0;
  n_b    (0)  = nb.GetY    (0);
  w_t    (0)  = wt.GetY    (0);
  w_p    (0)  = - wp.GetY    (0);
  w_E    (0)  = wE.GetY    (0);
  n_I    (0)  = nI.GetY    (0);
  dn_Idr (0)  = nI.GetdYdX (0) * fac0;
  T_I    (0)  = Ti.GetY    (0);
  dT_Idr (0)  = Ti.GetdYdX (0) * fac0;

  dn_edP1 (0) = ne.GetdYdX (0);
  dT_edP1 (0) = Te.GetdYdX (0);
  dn_idP1 (0) = ni.GetdYdX (0);
  dT_idP1 (0) = Ti.GetdYdX (0);
  
  n_e    (NPSI-1) = ne.GetY    (ne.GetN()-1);
  dn_edr (NPSI-1) = ne.GetdYdX (ne.GetN()-1) * fac1;
  T_e    (NPSI-1) = Te.GetY    (Te.GetN()-1);
  dT_edr (NPSI-1) = Te.GetdYdX (Te.GetN()-1) * fac1;
  n_i    (NPSI-1) = ni.GetY    (ni.GetN()-1);
  dn_idr (NPSI-1) = ni.GetdYdX (ni.GetN()-1) * fac1;
  T_i    (NPSI-1) = Ti.GetY    (Ti.GetN()-1);
  dT_idr (NPSI-1) = Ti.GetdYdX (Ti.GetN()-1) * fac1;
  n_b    (NPSI-1) = nb.GetY    (nb.GetN()-1);
  w_t    (NPSI-1) = wt.GetY    (wt.GetN()-1);
  w_p    (NPSI-1) =  - wp.GetY    (wp.GetN()-1);
  w_E    (NPSI-1) = wE.GetY    (wE.GetN()-1);
  n_I    (NPSI-1) = nI.GetY    (nI.GetN()-1);
  dn_Idr (NPSI-1) = nI.GetdYdX (nI.GetN()-1) * fac1;
  T_I    (NPSI-1) = Ti.GetY    (Ti.GetN()-1);
  dT_Idr (NPSI-1) = Ti.GetdYdX (Ti.GetN()-1) * fac1;

  dn_edP1 (NPSI-1) = ne.GetdYdX (ne.GetN()-1);
  dT_edP1 (NPSI-1) = Te.GetdYdX (Te.GetN()-1);
  dn_idP1 (NPSI-1) = ni.GetdYdX (ni.GetN()-1);
  dT_idP1 (NPSI-1) = Ti.GetdYdX (Ti.GetN()-1);
  
  for (int j = 1; j < NPSI-1; j++)
    {
      n_e    (j) = InterpolateField (ne, psi (j), 0);
      dn_edr (j) = InterpolateField (ne, psi (j), 1) * dpsidr (j) /a;
      T_e    (j) = InterpolateField (Te, psi (j), 0);
      dT_edr (j) = InterpolateField (Te, psi (j), 1) * dpsidr (j) /a;
      n_i    (j) = InterpolateField (ni, psi (j), 0);
      dn_idr (j) = InterpolateField (ni, psi (j), 1) * dpsidr (j) /a;
      T_i    (j) = InterpolateField (Ti, psi (j), 0);
      dT_idr (j) = InterpolateField (Ti, psi (j), 1) * dpsidr (j) /a;
      n_b    (j) = InterpolateField (nb, psi (j), 0);
      w_p    (j) = - InterpolateField (wp, psi (j), 0);
      w_t    (j) = InterpolateField (wt, psi (j), 0);
      w_E    (j) = InterpolateField (wE, psi (j), 0);
      n_I    (j) = InterpolateField (nI, psi (j), 0);
      dn_Idr (j) = InterpolateField (nI, psi (j), 1) * dpsidr (j) /a;
      T_I    (j) = InterpolateField (Ti, psi (j), 0);
      dT_Idr (j) = InterpolateField (Ti, psi (j), 1) * dpsidr (j) /a;

      dn_edP1 (j) = InterpolateField (ne, psi (j), 1);
      dT_edP1 (j) = InterpolateField (Te, psi (j), 1);
      dn_idP1 (j) = InterpolateField (ni, psi (j), 1);
      dT_idP1 (j) = InterpolateField (Ti, psi (j), 1);
     }

  // Calculate higher derivatives of density and temperature profiles
  for (int i = 0; i < NSMOOTH; i++)
    {
      Smoothing (NPSI, dn_edP1);
      Smoothing (NPSI, dT_edP1);
      Smoothing (NPSI, dn_idP1);
      Smoothing (NPSI, dT_idP1);
    }

  for (int j = 0; j < NPSI; j++)
    {
      dn_edP2 (j) = Interpolate (NPSI, psi, dn_edP1, psi (j), 1);
      dT_edP2 (j) = Interpolate (NPSI, psi, dT_edP1, psi (j), 1);
      dn_idP2 (j) = Interpolate (NPSI, psi, dn_idP1, psi (j), 1);
      dT_idP2 (j) = Interpolate (NPSI, psi, dT_idP1, psi (j), 1);
    }

  for (int i = 0; i < NSMOOTH; i++)
    {
      Smoothing (NPSI, dn_edP2);
      Smoothing (NPSI, dT_edP2);
      Smoothing (NPSI, dn_idP2);
      Smoothing (NPSI, dT_idP2);
    }
  
  for (int j = 0; j < NPSI; j++)
    {
      dn_edP3 (j) = Interpolate (NPSI, psi, dn_edP2, psi (j), 1);
      dT_edP3 (j) = Interpolate (NPSI, psi, dT_edP2, psi (j), 1);
      dn_idP3 (j) = Interpolate (NPSI, psi, dn_idP2, psi (j), 1);
      dT_idP3 (j) = Interpolate (NPSI, psi, dT_idP2, psi (j), 1);
    }

   for (int i = 0; i < NSMOOTH; i++)
    {
      Smoothing (NPSI, dn_edP3);
      Smoothing (NPSI, dT_edP3);
      Smoothing (NPSI, dn_idP3);
      Smoothing (NPSI, dT_idP3);
    }
  
  // Derived parameters
  for (int j = 0; j < NPSI; j++)
    {
      Quasi (j) = (ZI * (n_i (j) + n_b (j)) + ZII * n_I (j) - n_e (j)) /n_e (j);
      Z_eff (j) = (ZI*ZI * n_i (j) + ZII*ZII * n_I (j)) /n_e (j);
      alpha (j) = ZII*ZII * n_I (j) /ZI/ZI /n_i (j);
      if (NTYPE == 0)
	n_n (j) = NN * exp ((rr(j) - 1.) /(LN /a));
      else if (NTYPE == 1)
	n_n (j) = NN /(1. + (rr(j) - 1.) * (rr(j) - 1.) /(LN /a) /(LN /a));
    }
  
  // Interpolate cFiles
  if (INTC != 0)
    {
      // Save pwd
      char pwd[MAXFILENAMELENGTH];
      getcwd (pwd, MAXFILENAMELENGTH);

      // Remove cFile
      CallSystem ("rm -rf Inputs/cFile");

      // Get cFiles directory
      char cFileDir[MAXFILENAMELENGTH];
      CallSystem ("greadlink -f Inputs/cFiles > cFileDir");
      FILE* pfd = OpenFiler ("cFileDir");
      fscanf (pfd, "%s", cFileDir);
      fclose (pfd);
      CallSystem ("rm cFileDir");
  
      // Read cFile data
      char           Basename[MAXFILENAMELENGTH];
      char           Filename[MAXFILENAMELENGTH];
      char           filename[MAXFILENAMELENGTH];
      vector<string> cFileName;
      double         filetime;
      vector<double> cFileTime;
      int            cFileNumber;
  
      printf ("Reading cFile data:\n");

      chdir (cFileDir);
      getcwd (Basename, MAXFILENAMELENGTH);
      strcat (Basename, "/");
   
      FILE* file = OpenFiler ((char*) "Index");
   
      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  strcpy (Filename, Basename);
	  strcat (Filename, filename);
	  
	  cFileName.push_back (Filename);
	  cFileTime.push_back (filetime);
	}
      cFileNumber = cFileTime.size ();

      fclose (file);
      chdir (pwd);

      // Interpolate cFiles
      cFileInterp (cFileName, cFileTime, cFileNumber, TIME);
    }
  
  // Read perpendicular momentum diffusivity data from cFile
  printf ("Reading perpendicular momentum diffusivity data:\n");
  cFileRead ();

  // Interpolate perpendicular diffusivity onto existing psi grid
  chip.resize (NPSI);
  chie.resize (NPSI);
  chin.resize (NPSI);
  chii.resize (NPSI);
    
  chip (0)      = Chip.GetY (0);
  chie (0)      = Chie.GetY (0);
  chin (0)      = Chin.GetY (0);
  chii (0)      = Chii.GetY (0);
  chip (NPSI-1) = Chip.GetY (Chip.GetN()-1);
  chie (NPSI-1) = Chie.GetY (Chie.GetN()-1);
  chin (NPSI-1) = Chin.GetY (Chin.GetN()-1);
  chii (NPSI-1) = Chii.GetY (Chii.GetN()-1);

  for (int j = 1; j < NPSI-1; j++)
    {
      chip (j) = InterpolateField (Chip, psi (j), 0);
      chie (j) = InterpolateField (Chie, psi (j), 0);
      chin (j) = InterpolateField (Chin, psi (j), 0);
      chii (j) = InterpolateField (Chii, psi (j), 0);

      if (chip (j) < 0.)
	{
	  printf ("NEOCLASSICAL:: Warning - chip(%3d) < 0.\n", j);
	  chip (j) = 1.e-3;
	}
      if (chie (j) < 0.)
	{
	  printf ("NEOCLASSICAL:: Warning - chie(%3d) < 0.\n", j);
	  chie (j) = 1.e-3;
	}
      if (chin (j) < 0.)
	{
	  printf ("NEOCLASSICAL:: Warning - chin(%3d) < 0.\n", j);
	  chin (j) = 1.e-3; 
	}
      if (chii (j) < 0.)
	{
	  printf ("NEOCLASSICAL:: Warning - chii(%3d) < 0.\n", j);
	  chii (j) = 1.e-3; 
	}
    }
}

// #################################################
// Calculate derived quantities at rational surfaces
// #################################################
void Neoclassical::Get_Derived ()
{
  // -------------------------------------------
  // Calculate profile data at rational surfaces
  // -------------------------------------------
  rho0  = (AI * (n_i (0) + n_b (0)) + AII * n_I (0)) * m_p;
  tau_A = sqrt (mu_0 * rho0 * a*a /B_0/B_0);
  P0    = n_i (0) * T_i (0) + n_I (0) * T_I (0) + n_e (0) * T_e (0);
   
  nek.resize    (nres); dnedrk.resize (nres);
  Tek.resize    (nres); dTedrk.resize (nres);
  nik.resize    (nres); dnidrk.resize (nres);
  Tik.resize    (nres); dTidrk.resize (nres);
  nIk.resize    (nres); dnIdrk.resize (nres);
  TIk.resize    (nres); dTIdrk.resize (nres);
  wEk.resize    (nres); wtk.resize    (nres);
  nbk.resize    (nres); Zeffk.resize  (nres);
  alphak.resize (nres); rhok.resize   (nres);
  NNk.resize    (nres); chipk.resize  (nres);
  chiek.resize  (nres); chink.resize  (nres);
  chiik.resize  (nres); Zeffik.resize (nres);
  ZeffIk.resize (nres); wpk.resize    (nres);

  dnedP1k.resize (nres); dTedP1k.resize (nres); dnidP1k.resize (nres); dTidP1k.resize (nres);
  dnedP2k.resize (nres); dTedP2k.resize (nres); dnidP2k.resize (nres); dTidP2k.resize (nres);
  dnedP3k.resize (nres); dTedP3k.resize (nres); dnidP3k.resize (nres); dTidP3k.resize (nres);

  Factor1.resize  (nres); Factor2.resize  (nres); Factor3.resize  (nres);
  Factor4.resize  (nres); Factor5.resize  (nres); Factor6.resize  (nres);
  Factor7.resize  (nres); Factor8.resize  (nres); Factor9.resize  (nres);
  Factor10.resize (nres); Factor11.resize (nres); Factor12.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      nek    (j) = Interpolate (NPSI, rr, n_e,    rk (j), 0);
      dnedrk (j) = Interpolate (NPSI, rr, dn_edr, rk (j), 0);
      Tek    (j) = Interpolate (NPSI, rr, T_e,    rk (j), 0);
      dTedrk (j) = Interpolate (NPSI, rr, dT_edr, rk (j), 0);
      nik    (j) = Interpolate (NPSI, rr, n_i,    rk (j), 0);
      dnidrk (j) = Interpolate (NPSI, rr, dn_idr, rk (j), 0);
      Tik    (j) = Interpolate (NPSI, rr, T_i,    rk (j), 0);
      dTidrk (j) = Interpolate (NPSI, rr, dT_idr, rk (j), 0);
      nIk    (j) = Interpolate (NPSI, rr, n_I,    rk (j), 0);
      dnIdrk (j) = Interpolate (NPSI, rr, dn_Idr, rk (j), 0);
      TIk    (j) = Interpolate (NPSI, rr, T_I,    rk (j), 0);
      dTIdrk (j) = Interpolate (NPSI, rr, dT_Idr, rk (j), 0);
      wtk    (j) = Interpolate (NPSI, rr, w_t,    rk (j), 0);
      wpk    (j) = Interpolate (NPSI, rr, w_p,    rk (j), 0);
      wEk    (j) = Interpolate (NPSI, rr, w_E,    rk (j), 0);
      nbk    (j) = Interpolate (NPSI, rr, n_b,    rk (j), 0);
      Zeffk  (j) = Interpolate (NPSI, rr, Z_eff,  rk (j), 0);
      alphak (j) = Interpolate (NPSI, rr, alpha,  rk (j), 0);
      chipk  (j) = Interpolate (NPSI, rr, chip,   rk (j), 0);
      chiek  (j) = Interpolate (NPSI, rr, chie,   rk (j), 0);
      chink  (j) = Interpolate (NPSI, rr, chin,   rk (j), 0);
      chiik  (j) = Interpolate (NPSI, rr, chii,   rk (j), 0);

      if (chipk (j) < DMIN)
	chipk (j) = DMIN;
      if (chipk (j) > DMAX)
	chipk (j) = DMAX;
      if (chiek (j) < DMIN)
	chiek (j) = DMIN;
      if (chiek (j) > DMAX)
	chiek (j) = DMAX;
      if (chink (j) < DMIN)
	chink (j) = DMIN;
      if (chink (j) > DMAX)
	chink (j) = DMAX;
      if (chiik (j) < DMIN)
	chiik (j) = DMIN;
      if (chiik (j) > DMAX)
	chiik (j) = DMAX;

      rhok   (j) = (AI * (nik (j) + nbk (j)) + AII * nIk (j)) * m_p /rho0;
      Zeffik (j) =       (ZII - Zeffk (j)) /(ZII - 1.);
      ZeffIk (j) = ZII * (Zeffk (j) - 1.)  /(ZII - 1.);
      
      // Prevent negative diffusivities
      if (chipk (j) < 0.)
	chipk (j) = 1.;
      if (chiek (j) < 0.)
	chiek (j) = 1.;
      if (chiik (j) < 0.)
	chiik (j) = 1.;
      if (chink (j) < 0.)
	chink (j) = 1.;
  
      dnedP1k (j) = Interpolate (NPSI, psi, dn_edP1, rk (j), 0);
      dTedP1k (j) = Interpolate (NPSI, psi, dT_edP1, rk (j), 0);
      dnidP1k (j) = Interpolate (NPSI, psi, dn_idP1, rk (j), 0);
      dTidP1k (j) = Interpolate (NPSI, psi, dT_idP1, rk (j), 0);
      dnedP2k (j) = Interpolate (NPSI, psi, dn_edP2, rk (j), 0);
      dTedP2k (j) = Interpolate (NPSI, psi, dT_edP2, rk (j), 0);
      dnidP2k (j) = Interpolate (NPSI, psi, dn_idP2, rk (j), 0);
      dTidP2k (j) = Interpolate (NPSI, psi, dT_idP2, rk (j), 0);
      dnedP3k (j) = Interpolate (NPSI, psi, dn_edP3, rk (j), 0);
      dTedP3k (j) = Interpolate (NPSI, psi, dT_edP3, rk (j), 0);
      dnidP3k (j) = Interpolate (NPSI, psi, dn_idP3, rk (j), 0);
      dTidP3k (j) = Interpolate (NPSI, psi, dT_idP3, rk (j), 0);

      Factor1  (j) = dnedP1k (j) * Tek (j);
      Factor2  (j) = dTedP1k (j) * nek (j);
      Factor3  (j) = dnidP1k (j) * Tik (j);
      Factor4  (j) = dTidP1k (j) * nik (j);
      
      Factor5  (j) = - A2 (j) * (dnedP1k (j) * dTedP1k (j) + dnedP2k (j) * Tek (j)) /8.;
      Factor6  (j) = - A2 (j) * (dnedP1k (j) * dTedP1k (j) + dTedP2k (j) * nek (j)) /8.;
      Factor7  (j) = - A2 (j) * (dnidP1k (j) * dTidP1k (j) + dnidP2k (j) * Tik (j)) /8.;
      Factor8  (j) = - A2 (j) * (dnidP1k (j) * dTidP1k (j) + dTidP2k (j) * nik (j)) /8.;

      Factor9  (j) = dnedP3k (j) * Tek (j) /24.;
      Factor10 (j) = dTedP3k (j) * nek (j) /24.;
      Factor11 (j) = dnidP3k (j) * Tik (j) /24.;
      Factor12 (j) = dTidP3k (j) * nik (j) /24.;
      
      if (NTYPE == 0)
	NNk (j) = NN * exp ((rk(j) - 1.) /(LN /a));
      else if (NTYPE == 1)
	NNk (j) = NN / (1. + (rk(j) - 1.) * (rk(j) - 1.) /(LN /a) /(LN /a));
 
      printf ("m = %3d r = %10.3e ne = %10.3e Te = %10.3e ni = %10.3e Ti = %10.3e nI = %10.3e TI = %10.3e nb = %10.3e wt = %10.3e wp = %10.3e  wE = %10.3e Z_eff = %10.3e NN = %10.3e rho = %10.3e chip = %10.3e chie = %10.3e chin = %10.3e chii = %10.3e\n",
	      mk(j), rk(j), nek(j)/1.e19, Tek(j)/e/1.e3, nik(j)/1.e19, Tik(j)/e/1.e3, nIk(j)/1.e19, TIk(j)/e/1.e3, nbk(j)/1.e19, wtk(j)/1.e3, wpk(j)/1.e3, wEk(j)/1.e3, Zeffk(j), NNk(j)/1.e19, rhok(j), chipk(j), chiek(j), chink(j), chii(j));
    }

  if (IMPURITY == 0)
    for (int j = 0; j < nres; j++)
      {
	Zeffk  (j) = ZI;
	alphak (j) = 1.e-15;
      }

  if (NEUTRAL == 0)
    for (int j = 0; j < nres; j++)
      NNk (j) = 1.e-15;

  // -------------------------------------------------------------------------
  // Calculate thermal velocities and transit frequencies at rational surfaces
  // -------------------------------------------------------------------------
  v_T_ek.resize     (nres); v_T_ik.resize     (nres); v_T_Ik.resize     (nres);     
  omega_t_ek.resize (nres); omega_t_ik.resize (nres); omega_t_Ik.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      v_T_ek     (j) = sqrt (2. * Tek (j)      /m_e);
      v_T_ik     (j) = sqrt (2. * Tik (j) /AI  /m_p);
      v_T_Ik     (j) = sqrt (2. * TIk (j) /AII /m_p);
      omega_t_ek (j) = v_T_ek (j) * gmk (j) * Ktk (j);
      omega_t_ik (j) = v_T_ik (j) * gmk (j) * Ktk (j);
      omega_t_Ik (j) = v_T_Ik (j) * gmk (j) * Ktk (j);
    }

  // ----------------------------------------------------
  // Calculate collision frequencies at rational surfaces
  // ----------------------------------------------------
  nu_eek.resize (nres); nu_iik.resize (nres); nu_IIk.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      nu_eek (j) = (4./3./sqrt(M_PI)) * (4.*M_PI) * nek (j) * e*e*e*e * COULOMB
	/(4.*M_PI*epsilon_0)/(4.*M_PI*epsilon_0) /m_e/m_e /v_T_ek (j)/v_T_ek (j)/v_T_ek (j);
      nu_iik (j) = (4./3./sqrt(M_PI)) * (4.*M_PI) * nik (j) * e*e*e*e * ZI*ZI*ZI*ZI * COULOMB
	/(4.*M_PI*epsilon_0)/(4.*M_PI*epsilon_0) /(AI*m_p) /(AI*m_p) /v_T_ik (j)/v_T_ik (j)/v_T_ik (j);
      nu_IIk (j) = (4./3./sqrt(M_PI)) * (4.*M_PI) * nIk (j) * e*e*e*e * ZII*ZII*ZII*ZII * COULOMB
	/(4.*M_PI*epsilon_0)/(4.*M_PI*epsilon_0) /(AII*m_p) /(AII*m_p) /v_T_Ik (j)/v_T_Ik (j)/v_T_Ik (j);
    }

  // ----------------------------------------------------------------------------
  // Calculate critical island widths for profile flattening at rational surfaces
  // ----------------------------------------------------------------------------
  WcritTek.resize (nres); WcritTik.resize (nres); Wcritnek.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      double x               = 1.;
      double chi_para_e_brag = 1.581 /(1. + 0.2535 * Zeffk (j)) * v_T_ek (j) * v_T_ek (j) /nu_eek (j);
      for (int i = 0; i < 100; i++)
	{
	  double chi_para_e_max = 2. * R_0 * v_T_ek (j) /sqrt(M_PI) /double (ntor) /sk (j) /x;
	  double chi_para_e     = chi_para_e_brag * chi_para_e_max /(chi_para_e_brag + chi_para_e_max);

	  x = sqrt(8.) * pow (chiek (j) /chi_para_e, 0.25) * sqrt (R_0 /rk (j) /a /sk (j) /double (ntor));
	}
      WcritTek (j) = x * rk (j) * a;

      x = 1.;
      double chi_para_i_brag = sqrt(2.) * 1.953 * v_T_ik (j) * v_T_ik (j) /nu_iik (j) /Zeffk (j);
      for (int i = 0; i < 100; i++)
	{
	  double chi_para_i_max = 2. * R_0 * v_T_ik (j) /sqrt(M_PI) /double (ntor) /sk (j) /x;
	  double chi_para_i     = chi_para_i_brag * chi_para_i_max /(chi_para_i_brag + chi_para_i_max);

	  x = sqrt(8.) * pow (chiik (j) /chi_para_i, 0.25) * sqrt (R_0 /rk (j) /a /sk (j) /double (ntor));
	}
      WcritTik (j) = x * rk (j) * a;

      x = 1.;
      for (int i = 0; i < 100; i++)
	{
	  double chi_para_i_max = 2. * R_0 * v_T_ik (j) /sqrt(M_PI) /double (ntor) /sk (j) /x;
	  double chi_para_i     = chi_para_i_brag * chi_para_i_max /(chi_para_i_brag + chi_para_i_max);

	  x = sqrt(8.) * pow (chink (j) /chi_para_i, 0.25) * sqrt (R_0 /rk (j) /a /sk (j) /double (ntor));
	}
      Wcritnek (j) = x * rk (j) * a;

      printf ("m = %3d  WcritTe/a = %10.3e  WcritTi/a = %10.3e  Wcritne/a = %10.3e\n",
	      mk (j), WcritTek (j)/a, WcritTik (j)/a, Wcritnek (j)/a);
    }

  // ------------------------------------------------------
  // Calculate diamagnetic frequencies at rational surfaces
  // ------------------------------------------------------
  eta_ek.resize   (nres); eta_ik.resize   (nres); eta_Ik.resize   (nres);
  w_ast_ek.resize (nres); w_ast_ik.resize (nres); w_ast_Ik.resize (nres);
  w_betak.resize  (nres); w_Omegk.resize  (nres); 

  for (int j = 0; j < nres; j++)
    {
      eta_ek   (j) = nek (j) * dTedrk (j) /dnedrk (j) /Tek (j);
      eta_ik   (j) = nik (j) * dTidrk (j) /dnidrk (j) /Tik (j);		      
      eta_Ik   (j) = nIk (j) * dTIdrk (j) /dnIdrk (j) /TIk (j);
      w_ast_ek (j) =   (qk (j) /gk (j)) * Tek (j) * dnedrk (j) * (1. + eta_ek (j)) /e      /nek (j) /(a*rk (j)) /fabs (B_0);
      w_ast_ik (j) = - (qk (j) /gk (j)) * Tik (j) * dnidrk (j) * (1. + eta_ik (j)) /e /ZI  /nik (j) /(a*rk (j)) /fabs (B_0);
      w_ast_Ik (j) = - (qk (j) /gk (j)) * TIk (j) * dnIdrk (j) * (1. + eta_Ik (j)) /e /ZII /nIk (j) /(a*rk (j)) /fabs (B_0);
      w_betak  (j) = sk (j) * gk (j) * fabs (B_0) /mu_0 /nek (j) /R_0/R_0 /e /qk(j);
      w_Omegk  (j) = sk (j) * qk (j) * e * gk (j) * fabs (B_0) /AI /m_p;

      printf ("m = %3d wE = %10.3e w_ast_e = %10.3e w_ast_i = %10.3e w_ast_I = %10.3e w_beta = %10.3e w_Omega = %10.3e eta_e = %10.3e eta_i = %10.3e eta_I = %10.3e\n",
	      mk (j), wEk (j) /1.e3, w_ast_ek (j) /1.e3, w_ast_ik (j) /1.e3, w_ast_Ik (j) /1.e3, w_betak (j) /1.e3, w_Omegk (j) /1.e3, eta_ek (j), eta_ik(j), eta_Ik(j));
    }

  // -----------------------------------------
  // Calculate timescales at rational surfaces
  // -----------------------------------------
  rho_sk.resize (nres); tau_Ak.resize (nres); tau_Hk.resize (nres); tau_Rk.resize (nres); tau_Mk.resize (nres); tau_thk.resize (nres); tau_cxk.resize (nres);
 
  for (int j = 0; j < nres; j++)
    {
      rho_sk  (j) = sqrt ((5./3.) * AI * m_p * (Tek (j) + (nik (j) /nek (j)) * Tik (j) + (nIk (j) /nek (j)) * TIk (j))) /e /fabs (B_0) /gk (j);

      tau_Ak  (j) = tau_A * sqrt (rhok(j));
      tau_Hk  (j) = tau_A * sqrt (rhok(j)) * R_0 /a /gk(j) /sk(j) /double (ntor);
      tau_Rk  (j) = mu_0 * a*a * rk(j)*rk(j) * nek (j) * e*e /nu_eek (j) /m_e;
      tau_Mk  (j) = a*a /chipk(j);
      tau_thk (j) = 1. /nu_iik (j) / (1. + (qk(j) * R_0 /rk(j) /a) * (qk(j) * R_0 /rk(j) /a) /akk(j));
      tau_cxk (j) = 1. /NNk(j) /SVN;

      printf ("m = %3d  P0 = %10.3e  tau_A = %10.3e  tau_R = %10.3e  tau_M = %10.3e  tau_th = %10.3e  tau_cx = %10.3e\n",
	      mk (j), P0 /1.e19/e/1.e3, tau_A, tau_Rk (j), tau_Mk (j), tau_thk (j), tau_cxk (j));
    }

  // --------------------------------------------------------
  // Calculate collisionality parameters at rational surfaces
  // --------------------------------------------------------
  gt.resize      (nres);
  nu_P_e.resize  (nres); nu_P_i.resize  (nres); nu_P_I.resize  (nres);
  nu_PS_e.resize (nres); nu_PS_i.resize (nres); nu_PS_I.resize (nres);
  x_iI.resize    (nres); x_Ii.resize    (nres);

  for (int j = 0; j < nres; j++)
    {
      gt      (j) = (1. - fck (j)) /fck (j);
      nu_P_e  (j) = Kastk (j) * gt (j) * nu_eek (j) /omega_t_ek (j);
      nu_P_i  (j) = Kastk (j) * gt (j) * nu_iik (j) /omega_t_ik (j);
      nu_P_I  (j) = Kastk (j) * gt (j) * nu_IIk (j) /omega_t_Ik (j);
      nu_PS_e (j) = (5.*M_PI/8.) * nu_eek (j) /omega_t_ek (j);
      nu_PS_i (j) = (5.*M_PI/8.) * nu_iik (j) /omega_t_ik (j);
      nu_PS_I (j) = (5.*M_PI/8.) * nu_IIk (j) /omega_t_Ik (j);
      x_iI    (j) = v_T_Ik (j) /v_T_ik (j);
      x_Ii    (j) = v_T_ik (j) /v_T_Ik (j);

      printf ("m = %3d nu_P_e = %10.3e nu_P_i = %10.3e nu_P_I = %10.3e nu_PS_e = %10.3e nu_PS_i = %10.3e nu_PS_I = %10.3e\n",
	      mk (j), nu_P_e (j), nu_P_i (j), nu_P_I (j), nu_PS_e (j), nu_PS_i (j), nu_PS_I (j));
    }
}

// #######################################################
// Calculate neoclassical viscosities at rational surfaces
// #######################################################
void Neoclassical::Get_Viscosities ()
{
  mu_00_i.resize (nres); mu_01_i.resize (nres); mu_11_i.resize (nres);
  mu_00_I.resize (nres); mu_01_I.resize (nres); mu_11_I.resize (nres);
  mu_00_e.resize (nres); mu_01_e.resize (nres); mu_11_e.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      jj = j;

      double          x, h = h0, t_err;
      int             rept, step = 0, skip = 0; count = 0; flag = 0;
      Array<double,1> y (9);

      x = xmin;
      for (int i = 0; i < 9; i++)
	y(i) = 0.;

      do
	{
	  RK4Adaptive (x, y, h, t_err, acc, 2., rept, maxrept, hmin, hmax, 2, 0, NULL);
	}
      while (x < xmax);

      mu_00_i (j) = y (0);
      mu_01_i (j) = (5./2.) * y (0) - y (1);
      mu_11_i (j) = y (2) - 5. * y (1) + (25./4.) * y (0);

      mu_00_I (j) = y (3);
      mu_01_I (j) = (5./2.) * y (3) - y (4);
      mu_11_I (j) = y (5) - 5. * y (4) + (25./4.) * y (3);

      mu_00_e (j) = y (6);
      mu_01_e (j) = (5./2.) * y (6) - y (7);
      mu_11_e (j) = y (8) - 5. * y (7) + (25./4.) * y (6);

      printf ("m = %3d  mu_i = (%10.3e, %10.3e, %10.3e)  mu_I = (%10.3e, %10.3e, %10.3e)  mu_e = (%10.3e, %10.3e, %10.3e)\n",
	      mk (j), mu_00_i (j), mu_01_i (j), mu_11_i (j), mu_00_I (j), mu_01_I (j), mu_11_I (j), mu_00_e (j), mu_01_e (j), mu_11_e (j));
    }
}

// ######################################################
// Calculate neoclassical parameters at rational surfaces
// ######################################################
void Neoclassical::Get_Parameters ()
{
  L_ii_00.resize (nres); L_ii_01.resize (nres); L_iI_00.resize (nres);
  L_iI_01.resize (nres); L_Ii_00.resize (nres); L_Ii_01.resize (nres);
  L_II_00.resize (nres); L_II_01.resize (nres); G_ii_00.resize (nres);
  G_Ii_00.resize (nres); Q_00.resize    (nres); G_ei_00.resize (nres);
  L_ee_00.resize (nres); L_ee_01.resize (nres); L_ei_00.resize (nres);
  L_ei_01.resize (nres); L_eI_00.resize (nres); L_eI_01.resize (nres);
 
  // ........................
  // Calculate ion parameters
  // ........................
  for (int j = 0; j < nres; j++)
    {
      // [a] = [F^ii]
      double a00 =           alphak (j) * (1. + AI /AII) /pow (1. + x_iI (j)*x_iI (j), 1.5);
      double a01 = (3./2.) * alphak (j) * (1. + AI /AII) /pow (1. + x_iI (j)*x_iI (j), 2.5);
      double a10 = a01;
      double a11 =           alphak (j) * ((13./4.) + 4. * pow (x_iI (j), 2.) + (15./2.) * pow (x_iI (j), 4.))
	/pow (1. + x_iI (j)*x_iI (j), 2.5) + sqrt (2.);
      
      // [b] = [F^iI]
      double b00 = a00;
      double b01 = (3./2.)  * (Tik (j) /TIk (j)) * alphak (j) * (1. + AII /AI) /pow (1. + x_Ii (j)*x_Ii (j), 2.5) /x_iI (j);
      double b10 = a10;
      double b11 = (27./4.) * (Tik (j) /TIk (j)) * alphak (j) * pow (x_iI (j), 2.) /pow (1. + x_iI (j)*x_iI (j), 2.5);

      // [c] = [F^Ii]
      double c00 = a00;
      double c01 = a01;
      double c10 = b01;
      double c11 = (27./4.) * alphak (j) * pow (x_iI (j), 2.) /pow (1. + x_iI (j)*x_iI (j), 2.5);
											
      // [d] = [F^II]
      double d00 = a00;           
      double d01 = b01;
      double d10 = b01;
      double d11 = (Tik (j) /TIk (j)) * (alphak (j) * ((15./2.) + 4. * pow (x_iI (j), 2.) + (13./4.) * pow (x_iI (j), 4.))
					 /pow (1. + x_iI (j)*x_iI (j), 2.5) + sqrt (2.) * alphak (j)*alphak (j) * x_Ii (j));

      gsl_matrix* H    = gsl_matrix_alloc (4, 4);
      gsl_matrix* G    = gsl_matrix_alloc (4, 4);
      gsl_matrix* L    = gsl_matrix_alloc (4, 4);
      gsl_matrix* invH = gsl_matrix_alloc (4, 4);

      gsl_matrix* Lii  = gsl_matrix_alloc (2, 2);
      gsl_matrix* LiI  = gsl_matrix_alloc (2, 2);
      gsl_matrix* LIi  = gsl_matrix_alloc (2, 2);
      gsl_matrix* LII  = gsl_matrix_alloc (2, 2);

      gsl_matrix* Gii  = gsl_matrix_alloc (2, 2);
      gsl_matrix* GiI  = gsl_matrix_alloc (2, 2);
      gsl_matrix* GIi  = gsl_matrix_alloc (2, 2);
      gsl_matrix* GII  = gsl_matrix_alloc (2, 2);
      
      gsl_matrix_set (H, 0, 0, a00 + mu_00_i (j) + SVN * NNk (j)      /nu_iik (j) /YN);
      gsl_matrix_set (H, 0, 1, a01 + mu_01_i (j));
      gsl_matrix_set (H, 1, 0, a10 + mu_01_i (j));
      gsl_matrix_set (H, 1, 1, a11 + mu_11_i (j) + SVN * NNk (j) * EN /nu_iik (j) /YN);

      gsl_matrix_set (H, 0, 2, - b00);
      gsl_matrix_set (H, 0, 3, - b01);
      gsl_matrix_set (H, 1, 2, - b10);
      gsl_matrix_set (H, 1, 3, - b11);

      gsl_matrix_set (H, 2, 0, - c00);
      gsl_matrix_set (H, 2, 1, - c01);
      gsl_matrix_set (H, 3, 0, - c10);
      gsl_matrix_set (H, 3, 1, - c11);

      gsl_matrix_set (H, 2, 2, d00 + alphak (j)*alphak (j) * (Tik (j) / TIk (j)) * x_Ii (j) * mu_00_I (j));
      gsl_matrix_set (H, 2, 3, d01 + alphak (j)*alphak (j) * (Tik (j) / TIk (j)) * x_Ii (j) * mu_01_I (j));
      gsl_matrix_set (H, 3, 2, d10 + alphak (j)*alphak (j) * (Tik (j) / TIk (j)) * x_Ii (j) * mu_01_I (j));
      gsl_matrix_set (H, 3, 3, d11 + alphak (j)*alphak (j) * (Tik (j) / TIk (j)) * x_Ii (j) * mu_11_I (j));

      gsl_matrix_set (G, 0, 0, a00 + SVN * NNk (j)      /nu_iik (j));
      gsl_matrix_set (G, 0, 1, a01);
      gsl_matrix_set (G, 1, 0, a10);
      gsl_matrix_set (G, 1, 1, a11 + SVN * NNk (j) * EN /nu_iik (j));

      gsl_matrix_set (G, 0, 2, - b00);
      gsl_matrix_set (G, 0, 3, - b01);
      gsl_matrix_set (G, 1, 2, - b10);
      gsl_matrix_set (G, 1, 3, - b11);

      gsl_matrix_set (G, 2, 0, - c00);
      gsl_matrix_set (G, 2, 1, - c01);
      gsl_matrix_set (G, 3, 0, - c10);
      gsl_matrix_set (G, 3, 1, - c11);

      gsl_matrix_set (G, 2, 2, d00);
      gsl_matrix_set (G, 2, 3, d01);
      gsl_matrix_set (G, 3, 2, d10);
      gsl_matrix_set (G, 3, 3, d11);

      int s;

      gsl_permutation* p = gsl_permutation_alloc (4);

      gsl_linalg_LU_decomp (H, p, &s);

      gsl_linalg_LU_invert (H, p, invH);

      Matrix_Mult (invH, G, L, 4);

      L_ii_00 (j) = gsl_matrix_get (L, 0, 0);
      L_ii_01 (j) = gsl_matrix_get (L, 0, 1);
      L_iI_00 (j) = gsl_matrix_get (L, 0, 2);
      L_iI_01 (j) = gsl_matrix_get (L, 0, 3);
      L_Ii_00 (j) = gsl_matrix_get (L, 2, 0);
      L_Ii_01 (j) = gsl_matrix_get (L, 2, 1);
      L_II_00 (j) = gsl_matrix_get (L, 2, 2);
      L_II_01 (j) = gsl_matrix_get (L, 2, 3);
      G_ii_00 (j) = gsl_matrix_get (invH, 0, 0) * SVN * NNk (j) /nu_iik (j);
      G_Ii_00 (j) = gsl_matrix_get (invH, 2, 0) * SVN * NNk (j) /nu_iik (j);

      for (int is = 0; is < 2; is++)
	for (int js = 0; js < 2; js++)
	  {
	    gsl_matrix_set (Lii, is, js, gsl_matrix_get (L, is,     js));
	    gsl_matrix_set (LiI, is, js, gsl_matrix_get (L, is,     js + 2));
	    gsl_matrix_set (LIi, is, js, gsl_matrix_get (L, is + 2, js));
	    gsl_matrix_set (LII, is, js, gsl_matrix_get (L, is + 2, js + 2));

	    gsl_matrix_set (Gii, is, js, gsl_matrix_get (G, is,     js)     * SVN * NNk (j) /nu_iik (j));
	    gsl_matrix_set (GiI, is, js, gsl_matrix_get (G, is,     js + 2) * SVN * NNk (j) /nu_iik (j));
	    gsl_matrix_set (GIi, is, js, gsl_matrix_get (G, is + 2, js)     * SVN * NNk (j) /nu_iik (j));
	    gsl_matrix_set (GII, is, js, gsl_matrix_get (G, is + 2, js + 2) * SVN * NNk (j) /nu_iik (j));
	  }
            
      gsl_matrix_free (H); gsl_matrix_free (G); gsl_matrix_free (L); gsl_matrix_free (invH);
      gsl_permutation_free (p);

      // .............................
      // Calculate electron parameters
      // .............................
      gsl_matrix* Fee = gsl_matrix_alloc (2, 2);
      gsl_matrix* Fei = gsl_matrix_alloc (2, 2);
      gsl_matrix* FeI = gsl_matrix_alloc (2, 2);
      gsl_matrix* Gei = gsl_matrix_alloc (2, 2);
      gsl_matrix* Lee = gsl_matrix_alloc (2, 2);
      gsl_matrix* Lei = gsl_matrix_alloc (2, 2);
      gsl_matrix* LeI = gsl_matrix_alloc (2, 2);
      gsl_matrix* E   = gsl_matrix_alloc (2, 2);
      gsl_matrix* Q   = gsl_matrix_alloc (2, 2);

      gsl_matrix_set (Fee, 0, 0, Zeffk (j));
      gsl_matrix_set (Fee, 0, 1, (3./2.) * Zeffk (j));
      gsl_matrix_set (Fee, 1, 0, (3./2.) * Zeffk (j));
      gsl_matrix_set (Fee, 1, 1, sqrt(2.) + (13./4.) * Zeffk (j));

      gsl_matrix_set (Fei, 0, 0, Zeffik (j));
      gsl_matrix_set (Fei, 0, 1, 0.);
      gsl_matrix_set (Fei, 1, 0, (3./2.) * Zeffik (j));
      gsl_matrix_set (Fei, 1, 1, 0.);

      gsl_matrix_set (FeI, 0, 0, ZeffIk (j));
      gsl_matrix_set (FeI, 0, 1, 0.);
      gsl_matrix_set (FeI, 1, 0, (3./2.) * ZeffIk (j));
      gsl_matrix_set (FeI, 1, 1, 0.);

      gsl_matrix_set (E, 0, 0, Zeffk (j)                       + mu_00_e (j));
      gsl_matrix_set (E, 0, 1, (3./2.) * Zeffk (j)             + mu_01_e (j));
      gsl_matrix_set (E, 1, 0, (3./2.) * Zeffk (j)             + mu_01_e (j));
      gsl_matrix_set (E, 1, 1, sqrt(2.) + (13./4.) * Zeffk (j) + mu_11_e (j));

      gsl_permutation* pp = gsl_permutation_alloc (2);

      gsl_linalg_LU_decomp (E, pp, &s);

      gsl_linalg_LU_invert (E, pp, Q);

      Q_00 (j) = gsl_matrix_get (Q, 0, 0);
	    
      gsl_matrix* Temp1 = gsl_matrix_alloc (2, 2);
      gsl_matrix* Temp2 = gsl_matrix_alloc (2, 2);
      gsl_matrix* Temp3 = gsl_matrix_alloc (2, 2);

      Matrix_Mult (Fei,   Gii,   Temp1, 2);
      Matrix_Mult (FeI,   GIi,   Temp2, 2);
      Matrix_Add  (Temp1, Temp2, Temp3, 2);
      Matrix_Mult (Q,     Temp3, Gei,   2);

      Matrix_Mult (Q, Fee, Lee, 2);

      Matrix_Mult (Fei,   Lii,   Temp1, 2);
      Matrix_Sub  (Temp1, Fei,   Temp2, 2);
      Matrix_Mult (FeI,   LIi,   Temp3, 2);
      Matrix_Add  (Temp2, Temp3, Temp1, 2);
      Matrix_Mult (Q,     Temp1, Lei,   2);

      Matrix_Mult (FeI,   LII,   Temp1, 2);
      Matrix_Sub  (Temp1, FeI,   Temp2, 2);
      Matrix_Mult (Fei,   LiI,   Temp3, 2);
      Matrix_Add  (Temp2, Temp3, Temp1, 2);
      Matrix_Mult (Q,     Temp1, LeI,   2);

      G_ei_00 (j) = gsl_matrix_get (Gei, 0, 0);
      L_ee_00 (j) = gsl_matrix_get (Lee, 0, 0);
      L_ee_01 (j) = gsl_matrix_get (Lee, 0, 1);
      L_ei_00 (j) = gsl_matrix_get (Lei, 0, 0);
      L_ei_01 (j) = gsl_matrix_get (Lei, 0, 1);
      L_eI_00 (j) = gsl_matrix_get (LeI, 0, 0);
      L_eI_01 (j) = gsl_matrix_get (LeI, 0, 1);
      
      gsl_matrix_free (Temp1); gsl_matrix_free (Temp2); gsl_matrix_free (Temp3); 
      gsl_matrix_free (Lii);   gsl_matrix_free (LiI);   gsl_matrix_free (LIi);   gsl_matrix_free (LII);
      gsl_matrix_free (Gei);   gsl_matrix_free (Lee);   gsl_matrix_free (Lei);   gsl_matrix_free (LeI);
      gsl_matrix_free (Gii);   gsl_matrix_free (GiI);   gsl_matrix_free (GIi);   gsl_matrix_free (GII); 
      gsl_matrix_free (Fee);   gsl_matrix_free (Fei);   gsl_matrix_free (FeI);   gsl_matrix_free (E); gsl_matrix_free (Q);
      gsl_permutation_free (pp);
    }
  for (int j = 0; j < nres; j++)
    printf ("m = %3d  L_ii = (%10.3e, %10.3e)  L_iI = (%10.3e, %10.3e)  L_Ii = (%10.3e, %10.3e)  L_II = (%10.3e, %10.3e)  G_ii = %10.3e  G_Ii = %10.3e\n",
	    mk (j), L_ii_00 (j), L_ii_01 (j), L_iI_00 (j), L_iI_01 (j), L_Ii_00 (j), L_Ii_01 (j), L_II_00 (j), L_II_01 (j), G_ii_00 (j), G_Ii_00 (j));
  for (int j = 0; j < nres; j++)
    printf ("m = %3d  L_ee = (%10.3e, %10.3e)  L_ei = (%10.3e, %10.3e)  L_eI = (%10.3e, %10.3e)  G_ei = %10.3e  Q_00 = %10.3e\n",
	    mk (j), L_ee_00 (j), L_ee_01 (j), L_ei_00 (j), L_ei_01 (j), L_eI_00 (j), L_eI_01 (j), G_ei_00 (j), Q_00 (j));
  }

// #####################################################################
// Calculate natural island propagation frequencies at rational surfaces
// #####################################################################
void Neoclassical::Get_Frequencies ()
{
  w_linear.resize  (nres); w_nonlinear.resize (nres); w_EB.resize      (nres);
  w_E0k.resize     (nres); w_E1k.resize       (nres); w_E2k.resize     (nres);
  w_nc_I0k.resize  (nres); w_nc_I1k.resize    (nres); w_nc_I2k.resize  (nres); 
  w_pnc_I0k.resize (nres); w_pnc_I1k.resize   (nres); w_pnc_I2k.resize (nres);
  w_nc_ik.resize   (nres); w_nc_eek.resize    (nres); w_nc_eik.resize  (nres); 
  alpbek.resize    (nres); alpbik.resize      (nres); alpck.resize     (nres); 
  alppk.resize     (nres); rhothek.resize     (nres); rhothik.resize   (nres);

  // .....................
  // Calculate frequencies
  // .....................
  for (int j = 0; j < nres; j++)
    {
      // EXB = 0
      w_E0k     (j) = wEk (j);
      w_nc_I0k  (j) = - G_Ii_00 (j) * w_E0k (j) - L_II_00 (j) * w_ast_Ik (j) - L_Ii_00 (j) * w_ast_ik (j)
		      + L_II_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j)
	              + L_Ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j);
      w_pnc_I0k (j) = Kthek (j) * w_nc_I0k (j);

      // EXB = 1
      w_E1k     (j) = wtk (j) - w_ast_Ik (j) - wpk (j);
      w_nc_I1k  (j) = - G_Ii_00 (j) * w_E1k (j) - L_II_00 (j) * w_ast_Ik (j) - L_Ii_00 (j) * w_ast_ik (j)
		       + L_II_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j)
	               + L_Ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j);;
      w_pnc_I1k (j) = Kthek (j) * w_nc_I1k (j);
 
      // EXB = 2
      w_nc_I2k  (j)  = - L_II_00 (j) * w_ast_Ik (j) - L_Ii_00 (j) * w_ast_ik (j)
		       + L_II_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j)
	               + L_Ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j);
      w_E2k     (j)  = (wtk (j) - w_ast_Ik (j) - Kthek (j) * w_nc_I2k (j)) /(1. - Kthek (j) * G_Ii_00 (j));
      w_nc_I2k  (j) += - G_Ii_00 (j) * w_E2k (j);
      w_pnc_I2k (j)  = Kthek (j) * w_nc_I2k (j);

      if (EXB == 0)
	{
	  w_nc_ik (j)  = - G_ii_00 (j) * w_E0k (j) - L_ii_00 (j) * w_ast_ik (j) - L_iI_00 (j) * w_ast_Ik (j)
	                 + L_ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                 + L_iI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  w_nc_eek (j) = - L_ee_00 (j) * w_ast_ek (j)
	                 + L_ee_01 (j) * (eta_ek (j) /(1. + eta_ek (j))) * w_ast_ek (j);
	  w_nc_eik (j) = - G_ei_00 (j) * w_E0k (j) - L_ei_00 (j) * w_ast_ik (j) - L_eI_00 (j) * w_ast_Ik (j)
	                 + L_ei_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                 + L_eI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  
	  w_linear    (j) = - double (ntor) * (w_E0k (j) + w_ast_ek (j));
	  w_nonlinear (j) = - double (ntor) * (w_E0k (j) + w_ast_ik (j) + w_nc_ik (j));
	  w_EB        (j) = - double (ntor) * (w_E0k (j));
	}
      else if (EXB == 1)
	{
	  w_nc_ik (j)  = - G_ii_00 (j) * w_E1k (j) - L_ii_00 (j) * w_ast_ik (j) - L_iI_00 (j) * w_ast_Ik (j)
	                 + L_ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                 + L_iI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  w_nc_eek (j) = - L_ee_00 (j) * w_ast_ek (j)
	                 + L_ee_01 (j) * (eta_ek (j) /(1. + eta_ek (j))) * w_ast_ek (j);
	  w_nc_eik (j) = - G_ei_00 (j) * w_E1k (j) - L_ei_00 (j) * w_ast_ik (j) - L_eI_00 (j) * w_ast_Ik (j)
	                 + L_ei_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                 + L_eI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  
	  w_linear    (j) = - double (ntor) * (w_E1k (j) + w_ast_ek (j));
	  w_nonlinear (j) = - double (ntor) * (w_E1k (j) + w_ast_ik (j) + w_nc_ik (j));
	  w_EB        (j) = - double (ntor) * (w_E1k (j));
	}
      else if (EXB == 2)
	{
	  w_nc_ik (j)  = - G_ii_00 (j) * w_E2k (j) - L_ii_00 (j) * w_ast_ik (j) - L_iI_00 (j) * w_ast_Ik (j)
	                 + L_ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                 + L_iI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  w_nc_eek (j) = - L_ee_00 (j) * w_ast_ek (j)
	                 + L_ee_01 (j) * (eta_ek (j) /(1. + eta_ek (j))) * w_ast_ek (j);
	  w_nc_eik (j) = - G_ei_00 (j) * w_E2k (j) - L_ei_00 (j) * w_ast_ik (j) - L_eI_00 (j) * w_ast_Ik (j)
	                 + L_ei_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                 + L_eI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  
	  w_linear    (j) = - double (ntor) * (w_E2k (j) + w_ast_ek (j));
	  w_nonlinear (j) = - double (ntor) * (w_E2k (j) + w_ast_ik (j) + w_nc_ik (j));
	  w_EB        (j) = - double (ntor) * (w_E2k (j));
	}

      printf ("m = %3d  w_linear = %10.3e  w_nonlinear = %10.3e  w_EB = %10.3e  w_ast_ik = %10.3e  w_ast_ek = %10.3e  w_nc_i = %10.3e  w_nc_e = %10.3e\n",
	      mk (j), w_linear (j) /1.e3, w_nonlinear (j) /1.e3, w_EB (j) /1.e3, w_ast_ik (j) /1.e3, w_ast_ek (j) /1.e3, w_nc_ik (j) /1.e3, (w_nc_eek (j) + w_nc_eik (j)) /1.e3);
    }
  
  // ...........................................................
  // Calculate bootstrap, curvature, and polarization parameters
  // ...........................................................
  for (int j = 0; j < nres; j++)
    {
      alpbek(j) = - 2.*0.8227*1.58 * (w_ast_ek (j) + w_nc_eek (j)) /w_betak (j);
      if (EXB == 0)
	alpbik (j) = 2.*0.8227*1.58 * (- w_nc_eik (j) + (ZI * nik (j) /nek (j)) * (w_ast_ik (j) + w_nc_ik (j)) + (ZII * nIk (j) /nek (j)) * (w_ast_Ik (j) + w_nc_I0k (j))) /w_betak (j);
      else if (EXB == 1)
	alpbik (j) = 2.*0.8227*1.58 * (- w_nc_eik (j) + (ZI * nik (j) /nek (j)) * (w_ast_ik (j) + w_nc_ik (j)) + (ZII * nIk (j) /nek (j)) * (w_ast_Ik (j) + w_nc_I1k (j))) /w_betak (j);
      else  if (EXB == 2)
	alpbik (j) = 2.*0.8227*1.58 * (- w_nc_eik (j) + (ZI * nik (j) /nek (j)) * (w_ast_ik (j) + w_nc_ik (j)) + (ZII * nIk (j) /nek (j)) * (w_ast_Ik (j) + w_nc_I2k (j))) /w_betak (j);
      alpck   (j) =   2.*0.8227*1.58 * DR (j);
      alppk   (j) =   8.*0.8227*0.8227*0.8227*1.38 * (w_ast_ik (j) + w_nc_ik (j)) * w_nc_ik (j) /w_betak (j) /w_Omegk (j);
      rhothek (j) =   v_T_ek (j)      * m_e * qk (j) * R_0 /e /fabs(B_0) /gk (j) /rk (j) /a;
      rhothik (j) =   v_T_ik (j) * AI * m_p * qk (j) * R_0 /e /fabs(B_0) /gk (j) /rk (j) /a;

      printf ("m = %3d  alpha_b_e = %10.3e  alpha_b_i = %10.3e  alpha_c = %10.3e  alpha_p = %10.3e  rho_th_e/a = %10.3e  rho_th_i/a = %10.3e\n",
	      mk (j), alpbek (j), alpbik (j), alpck (j), alppk (j), rhothek (j)/a, rhothik (j)/a);
    }
}

// ##################################################
// Calculate linear layer widths at rational surfaces
// ##################################################
void Neoclassical::Get_LayerWidths ()
{
  // ---------------------------------
  // Calculate linear layer parameters
  // ---------------------------------
  Sk.resize  (nres); tauk.resize (nres); PEk.resize (nres); PMk.resize  (nres); Dk.resize (nres);
  QEk.resize (nres); Qek.resize  (nres); Qik.resize (nres); delk.resize (nres);
  
  for (int j = 0; j < nres; j++)
    {
      double tauE = a*a * rk (j)*rk (j) /chink (j);
      double tauM = a*a * rk (j)*rk (j) /chipk (j);

      Sk   (j) = Q_00 (j) * tau_Rk (j) /tau_Hk (j);
      tauk (j) = - w_ast_ek (j) /w_ast_ik (j);
      PEk  (j) = tau_Rk(j) /tauE;
      PMk  (j) = tau_Rk(j) /tauM;
      Dk   (j) = pow (Sk (j), 1./3.) * sqrt (tauk (j)/ (1. + tauk (j))) * rho_sk (j) /a/rk(j);
      if (EXB == 0)
	QEk (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_E0k (j) * tau_Hk (j);
      else if (EXB == 1)
	QEk (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_E1k (j) * tau_Hk (j);
      else if (EXB == 2)
	QEk (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_E2k (j) * tau_Hk (j);
      Qek  (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_ast_ek (j) * tau_Hk (j);
      Qik  (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_ast_ik (j) * tau_Hk (j);

      // Prevent 1 + tauk from taking negative value
      if (1. /tauk (j) < TAUMIN)
	tauk (j) = 1. /TAUMIN;

      printf ("m = %3d  S = %11.4e  tau = %10.3e  P_E = %10.3e  P_M = %10.3e  D = %10.3e  Q_E = %10.3e  Q_e = %10.3e  Q_i = %10.3e\n",
	      mk (j), Sk (j), tauk (j), PEk (j), PMk (j), Dk (j), QEk (j), Qek (j), Qik (j));
    }

  // -----------------------------
  // Calculate linear layer widths
  // -----------------------------
  for (int j = 0; j < nres; j++)
    {
      jres = j;
      
      double          x, h = h0, t_err, max = 1.e20;
      int             rept, step = 0, skip = 0; count = 0; flag = 1;
      Array<double,1> y (8);
      double          C    = PEk (j) / (1. + 1./tauk (j)) /Dk (j)/Dk (j);
      double          xmax = sqrt (2. * log (max) /sqrt (C));

      if (C < 0.)
	{
	  printf ("NEOCLASSICAL::Get_LayerWidths: Error - C < 0.\n");
	  exit (1);
	}

      x     = 0.;
      y (0) = 1.;
      y (1) = 0.;
      y (2) = 0.;
      y (3) = 0.;
      y (4) = 0.;
      y (5) = 0.;
      y (6) = 1.;
      y (7) = 0.;
      
      do
	{
	  RK4Adaptive (x, y, h, t_err, acc, 2., rept, maxrept, hmin, hmax, 2, 0, NULL);
	}
      while (x < xmax);

      gsl_complex Y1 = gsl_complex_rect (y (0), y (1));
      Y1             = gsl_complex_mul_real (Y1, sqrt (xmax) /exp (sqrt (C) * xmax*xmax /2.));
      gsl_complex Y2 = gsl_complex_rect (y (4), y (5));
      Y2             = gsl_complex_mul_real (Y2, sqrt (xmax) /exp (sqrt (C) * xmax*xmax /2.));

      gsl_complex c  = gsl_complex_div (Y1, Y2);

      delk (j) = M_PI * gsl_complex_abs (c) * a*rk(j) /pow (Sk (j), 1./3.);

      if (isnan (delk (j)))
	{
	  printf ("Neoclassical::Get_LayerWidths: Error - delk (%2d) = NaN\n", j);
	  exit (1);
	}
    }
}

// ####################################
// Function to write Stage3 NETCDF file
// ####################################
void Neoclassical::WriteStage3Netcdfc ()
{
  // Convert data from blitz++ array to c array
  double* psi_x         = new double[NPSI];
  double* n_e_x         = new double[NPSI];
  double* T_e_x         = new double[NPSI];
  double* n_i_x         = new double[NPSI];
  double* T_i_x         = new double[NPSI];
  double* n_I_x         = new double[NPSI];
  double* T_I_x         = new double[NPSI];
  double* w_t_x         = new double[NPSI];
  double* w_p_x         = new double[NPSI];
  double* w_E_x         = new double[NPSI];
  double* Z_eff_x       = new double[NPSI];
  double* n_n_x         = new double[NPSI];
  double* chip_x        = new double[NPSI];
  double* chie_x        = new double[NPSI];
  double* chin_x        = new double[NPSI];
  double* chii_x        = new double[NPSI];
  for (int i = 0; i < NPSI; i++)
    {
      psi_x[i]   = psi(i);
      n_e_x[i]   = n_e(i) /1.e19; 
      T_e_x[i]   = T_e(i) /e/1.e3;
      n_i_x[i]   = n_i(i) /1.e19; 
      T_i_x[i]   = T_i(i) /e/1.e3; 
      n_I_x[i]   = n_I(i) /1.e19; 
      T_I_x[i]   = T_I(i) /e/1.e3; 
      w_t_x[i]   = w_t(i) /1.e3;
      w_p_x[i]   = w_p(i) /1.e3;
      w_E_x[i]   = w_E(i) /1.e3;
      Z_eff_x[i] = Z_eff(i);
      n_n_x[i]   = n_n(i) /1.e19;
      chip_x[i]  = chip(i);
      chie_x[i]  = chie(i);
      chin_x[i]  = chin(i);
      chii_x[i]  = chii(i);
    }

  int*    mk_x          = new int[nres];
  double* PsiNk_x       = new double[nres];
  double* tau_Ak_x      = new double[nres];
  double* tau_Hk_x      = new double[nres];
  double* tau_Rk_x      = new double[nres];
  double* tau_Mk_x      = new double[nres];
  double* tau_thk_x     = new double[nres];
  double* tau_cxk_x     = new double[nres];
  double* w_linear_x    = new double[nres];
  double* w_nonlinear_x = new double[nres];
  double* w_EB_x        = new double[nres];
  double* rho_sk_x      = new double[nres];
  double* delk_x        = new double[nres];
  double* rhothek_x     = new double[nres];
  double* rhothik_x     = new double[nres];
  double* WcritTek_x    = new double[nres];
  double* WcritTik_x    = new double[nres];
  double* Wcritnek_x    = new double[nres];
  double* w_ast_Ik_x    = new double[nres];
  double* w_ast_ek_x    = new double[nres];
  double* w_ast_ik_x    = new double[nres];
  double* Sk_x          = new double[nres];
  double* tauk_x        = new double[nres];
  double* PEk_x         = new double[nres];
  double* PMk_x         = new double[nres];
  double* Dk_x          = new double[nres];
  double* QEk_x         = new double[nres];
  double* Qek_x         = new double[nres];
  double* Qik_x         = new double[nres];
  double* alpbek_x      = new double[nres];
  double* alpbik_x      = new double[nres];
  double* alpck_x       = new double[nres];
  double* alppk_x       = new double[nres];
  double* fac1_x        = new double[nres];
  double* fac2_x        = new double[nres];
  double* fac3_x        = new double[nres];
  double* fac4_x        = new double[nres];
  double* fac5_x        = new double[nres];
  double* fac6_x        = new double[nres];
  double* fac7_x        = new double[nres];
  double* fac8_x        = new double[nres];
  double* fac9_x        = new double[nres];
  double* fac10_x       = new double[nres];
  double* fac11_x       = new double[nres];
  double* fac12_x       = new double[nres];	
  double* w_tk_x        = new double[nres];
  double* w_pk_x        = new double[nres];
  double* w_EB0_x       = new double[nres];
  double* w_EB1_x       = new double[nres];
  double* w_EB2_x       = new double[nres];
  double* w_nc_I0_x     = new double[nres];
  double* w_nc_I1_x     = new double[nres];
  double* w_nc_I2_x     = new double[nres];
  double* w_pnc_I0_x    = new double[nres];
  double* w_pnc_I1_x    = new double[nres];
  double* w_pnc_I2_x    = new double[nres];
  double* w_nc_i_x      = new double[nres];
  double* w_nc_e_x      = new double[nres];
  double* w_th_i_x      = new double[nres];
  double* w_para_x      = new double[nres];  
  for (int i = 0; i < nres; i++)
    {
      mk_x[i]          = mk(i);
      PsiNk_x[i]       = PsiNk(i);
      tau_Ak_x[i]      = tau_Ak(i);
      tau_Hk_x[i]      = tau_Hk(i);
      tau_Rk_x[i]      = tau_Rk(i) * Q_00(i);
      tau_Mk_x[i]      = tau_Mk(i);
      tau_thk_x[i]     = tau_thk(i) /mu_00_i(i);
      tau_cxk_x[i]     = tau_cxk(i);
      w_linear_x[i]    = w_linear(i)    /1.e3;
      w_nonlinear_x[i] = w_nonlinear(i) /1.e3;
      w_EB_x[i]        = w_EB(i)        /1.e3;
      rho_sk_x[i]      = rho_sk(i)   /1.e-2;
      delk_x[i]        = delk(i)     /1.e-2;
      rhothek_x[i]     = rhothek(i)  /1.e-2;
      rhothik_x[i]     = rhothik(i)  /1.e-2;
      WcritTek_x[i]    = WcritTek(i) /1.e-2;
      WcritTik_x[i]    = WcritTik(i) /1.e-2;
      Wcritnek_x[i]    = Wcritnek(i) /1.e-2;
      w_tk_x[i]        = wtk(i) /1.e3;
      w_pk_x[i]        = wpk(i) /1.e3;
      w_EB0_x[i]       = w_E0k(i)    /1.e3;
      w_EB1_x[i]       = w_E1k(i)    /1.e3;
      w_EB2_x[i]       = w_E2k(i)    /1.e3;
      w_nc_I0_x[i]     = w_nc_I0k(i) /1.e3;
      w_nc_I1_x[i]     = w_nc_I1k(i) /1.e3;
      w_nc_I2_x[i]     = w_nc_I2k(i) /1.e3;  
      w_pnc_I0_x[i]    = w_pnc_I0k(i) /1.e3;
      w_pnc_I1_x[i]    = w_pnc_I1k(i) /1.e3;
      w_pnc_I2_x[i]    = w_pnc_I2k(i) /1.e3;
      w_ast_Ik_x[i]    = w_ast_Ik(i) /1.e3;
      w_ast_ek_x[i]    = w_ast_ek(i) /1.e3;
      w_ast_ik_x[i]    = w_ast_ik(i) /1.e3;
      Sk_x[i]          = Sk(i);
      tauk_x[i]        = tauk(i);
      PEk_x[i]         = PEk(i);
      PMk_x[i]         = PMk(i);
      Dk_x[i]          = Dk(i);
      QEk_x[i]         = QEk(i);
      Qek_x[i]         = Qek(i);
      Qik_x[i]         = Qik(i);
      alpbek_x[i]      = alpbek(i);
      alpbik_x[i]      = alpbik(i);
      alpck_x[i]       = alpck(i);
      alppk_x[i]       = alppk(i);
      fac1_x[i]        = Factor1(i)  /1.e19/e/1.e3;
      fac2_x[i]        = Factor2(i)  /1.e19/e/1.e3;
      fac3_x[i]        = Factor3(i)  /1.e19/e/1.e3;
      fac4_x[i]        = Factor4(i)  /1.e19/e/1.e3;
      fac5_x[i]        = Factor5(i)  /1.e19/e/1.e3;
      fac6_x[i]        = Factor6(i)  /1.e19/e/1.e3;
      fac7_x[i]        = Factor7(i)  /1.e19/e/1.e3;
      fac8_x[i]        = Factor8(i)  /1.e19/e/1.e3;
      fac9_x[i]        = Factor9(i)  /1.e19/e/1.e3;
      fac10_x[i]       = Factor10(i) /1.e19/e/1.e3;
      fac11_x[i]       = Factor11(i) /1.e19/e/1.e3;
      fac12_x[i]       = Factor12(i) /1.e19/e/1.e3;
      w_nc_i_x[i]      = w_nc_ik(i) /1.e3;
      w_nc_e_x[i]      = (w_nc_eek(i) + w_nc_eik(i)) /1.e3;
      w_th_i_x[i]      = Kthek(i) * w_nc_ik(i) /1.e3;
      if (EXB == 0)
	w_para_x[i]    = (- (w_ast_ek(i) + w_nc_eek(i) + w_nc_eik(i)) 
			  + (ZI  * nik(i) /nek(i)) * (w_ast_ik(i) + w_nc_ik(i))
			  + (ZII * nIk(i) /nek(i)) * (w_ast_Ik(i) + w_nc_I0k(i))) /1.e3;
      else if (EXB == 1)
	w_para_x[i]    = (- (w_ast_ek(i) + w_nc_eek(i) + w_nc_eik(i)) 
			  + (ZI  * nik(i) /nek(i)) * (w_ast_ik(i) + w_nc_ik(i))
			  + (ZII * nIk(i) /nek(i)) * (w_ast_Ik(i) + w_nc_I1k(i))) /1.e3;
      else if (EXB == 2)
	w_para_x[i]    = (- (w_ast_ek(i) + w_nc_eek(i) + w_nc_eik(i)) 
			  + (ZI  * nik(i) /nek(i)) * (w_ast_ik(i) + w_nc_ik(i))
			  + (ZII * nIk(i) /nek(i)) * (w_ast_Ik(i) + w_nc_I2k(i))) /1.e3;
    }

  int NISLAND = 200;
  Array<double,2> RHS1  (nres, NISLAND);
  Array<double,2> RHS2  (nres, NISLAND);
  Array<double,2> RHS3  (nres, NISLAND);
  Array<double,2> RHS1a (nres, NISLAND);
  Array<double,2> RHS2a (nres, NISLAND);
  Array<double,2> RHS3a (nres, NISLAND);

  for (int j = 0; j < nres; j++)
    for (int i = 0; i < NISLAND; i++)
      {
	double wpsi   = double (i+1) * WPSIMAX /double (NISLAND);
	double wr     = R_0 * wpsi /dPsidr (j);
	double wrh    = 0.8227 * wr           /2. /rk(j) /a;
	double wrce   = 0.8227 * WcritTek (j) /2. /rk(j) /a;
	double wrci   = 0.8227 * WcritTik (j) /2. /rk(j) /a;
	double wrcn   = 0.8227 * Wcritnek (j) /2. /rk(j) /a;
	double rhoe   = 0.8227 * rhothek  (j) /2. /rk(j) /a;
	double rhoi   = 0.8227 * rhothik  (j) /2. /rk(j) /a;
	double facbTe = (eta_ek (j) /(1. + eta_ek (j))) * wr /(wrce*wrce + rhoe*rhoe + wr*wr);
	double facbne = (1.         /(1. + eta_ek (j))) * wr /(wrcn*wrcn + rhoe*rhoe + wr*wr);
	double facbTi = (eta_ik (j) /(1. + eta_ik (j))) * wr /(wrci*wrci + rhoi*rhoi + wr*wr);
	double facbni = (1.         /(1. + eta_ik (j))) * wr /(wrcn*wrcn + rhoi*rhoi + wr*wr);
	double faccTe = (eta_ek (j) /(1. + eta_ek (j))) * wr /(wrce*wrce + wr*wr);
	double faccne = (1.         /(1. + eta_ek (j))) * wr /(wrcn*wrcn + wr*wr);
	double faccTi = (eta_ik (j) /(1. + eta_ik (j))) * wr /(wrci*wrci + wr*wr);
	double faccni = (1.         /(1. + eta_ik (j))) * wr /(wrcn*wrcn + wr*wr);
	double facc   =   (nek (j) /(nek (j) + nik (j))) * (faccTe + faccne)
	                + (nik (j) /(nek (j) + nik (j))) * (faccTi + faccni);
	double facpTi = (eta_ik (j) /(1. + eta_ik (j))) * wr /(wrci*wrci + wr*wr) /(wrci*wrci + wr*wr);
	double facpni = (1.         /(1. + eta_ik (j))) * wr /(wrcn*wrcn + wr*wr) /(wrcn*wrcn + wr*wr);
 
	double boot  =   alpbek (j) * (facbTe + facbne) + alpbik (j) * (facbTi + facbni);
	double curv  =   alpck  (j) * facc;
	double polz  =   alppk  (j) * (facpTi + facpni);
	double wall  =   Sigmaw (j) * Sigmaw (j) /Deltaw (j);
	double poem1 =   Poem1  (j) * (wr /rk (j) /a) * log (wr /rk(j) /a);
	double poem2 = - Poem2  (j) * wr /rk(j) /a;
	double poem3 = - Poem3  (j) * wr /rk(j) /a;
	double poem  = ((poem1 + poem2) *  wrce*wrce + poem3 * wr*wr) /(wrce*wrce + wr*wr);

	double rhs1 = - EEh (j, j) + poem;
	double rhs2 = rhs1 + boot + curv + polz;
	double rhs3 = rhs1 + boot + curv + polz + wall;

	double rhs1a = - EEh (j, j) + poem3;
	double rhs2a = rhs1a + boot + curv + polz;
	double rhs3a = rhs1a + boot + curv + polz + wall;

	RHS1 (j, i) = rhs1;
	RHS2 (j, i) = rhs2;
	RHS3 (j, i) = rhs3;

	RHS1a (j, i) = rhs1a;
	RHS2a (j, i) = rhs2a;
	RHS3a (j, i) = rhs3a;
      }

  double* w_x     = new double [NISLAND];
  double* rhs1_x  = new double [nres*NISLAND];
  double* rhs2_x  = new double [nres*NISLAND];
  double* rhs3_x  = new double [nres*NISLAND];
  double* rhs1a_x = new double [nres*NISLAND];
  double* rhs2a_x = new double [nres*NISLAND];
  double* rhs3a_x = new double [nres*NISLAND];
  for (int j = 0; j < nres; j++)
    for (int i = 0; i < NISLAND; i++)
      {
	rhs1_x [i + j*NISLAND] = RHS1  (j, i);
	rhs2_x [i + j*NISLAND] = RHS2  (j, i);
	rhs3_x [i + j*NISLAND] = RHS3  (j, i);
	rhs1a_x[i + j*NISLAND] = RHS1a (j, i);
	rhs2a_x[i + j*NISLAND] = RHS2a (j, i);
	rhs3a_x[i + j*NISLAND] = RHS3a (j, i);
      }

  for (int i = 0; i < NISLAND; i++)
    w_x[i] = double (i+1) * WPSIMAX /double (NISLAND);
  
  // Open file
  int err = 0, dataFile;
  err = nc_create ("Outputs/Stage3.nc", NC_CLOBBER, &dataFile);
 
  if (err != 0)
    {
      printf ("NEOCLASSICAL::WriteStage3Netcdfc: Error opening Outputs/Stage3.nc\n");
      exit (1);
    }

  // PsiN
  int PsiN_d, PsiN;
  err += nc_def_dim (dataFile, "N_psi", NPSI, &PsiN_d);
  err += nc_def_var (dataFile, "PsiN", NC_DOUBLE, 1, &PsiN_d, &PsiN);

  // n_e
  int n_e_y;
  err += nc_def_var (dataFile, "n_e", NC_DOUBLE, 1, &PsiN_d, &n_e_y);

  // T_e
  int T_e_y;
  err += nc_def_var (dataFile, "T_e", NC_DOUBLE, 1, &PsiN_d, &T_e_y);
  
  // n_i
  int n_i_y;
  err += nc_def_var (dataFile, "n_i", NC_DOUBLE, 1, &PsiN_d, &n_i_y);

  // T_i
  int T_i_y;
  err += nc_def_var (dataFile, "T_i", NC_DOUBLE, 1, &PsiN_d, &T_i_y);
  
  // n_I
  int n_I_y;
  err += nc_def_var (dataFile, "n_I", NC_DOUBLE, 1, &PsiN_d, &n_I_y);

  // T_I
  int T_I_y;
  err += nc_def_var (dataFile, "T_I", NC_DOUBLE, 1, &PsiN_d, &T_I_y);
  
  // Z_eff
  int Z_eff_y;
  err += nc_def_var (dataFile, "Z_eff", NC_DOUBLE, 1, &PsiN_d, &Z_eff_y);

  // n_n
  int n_n_y;
  err += nc_def_var (dataFile, "n_n", NC_DOUBLE, 1, &PsiN_d, &n_n_y);
  
  // w_t
  int w_t_y;
  err += nc_def_var (dataFile, "w_t", NC_DOUBLE, 1, &PsiN_d, &w_t_y);

  // w_p
  int w_p_y;
  err += nc_def_var (dataFile, "w_p", NC_DOUBLE, 1, &PsiN_d, &w_p_y);

  // w_E
  int w_E_y;
  err += nc_def_var (dataFile, "w_E", NC_DOUBLE, 1, &PsiN_d, &w_E_y);

  // chip
  int chip_y;
  err += nc_def_var (dataFile, "chi_p", NC_DOUBLE, 1, &PsiN_d, &chip_y);

  // chie
  int chie_y;
  err += nc_def_var (dataFile, "chi_e", NC_DOUBLE, 1, &PsiN_d, &chie_y);

  // chin
  int chin_y;
  err += nc_def_var (dataFile, "chi_n", NC_DOUBLE, 1, &PsiN_d, &chin_y);

  // chii
  int chii_y;
  err += nc_def_var (dataFile, "chi_i", NC_DOUBLE, 1, &PsiN_d, &chii_y);

  // mk
  int nres_d, mk_y;
  err += nc_def_dim (dataFile, "N_res", nres, &nres_d);
  err += nc_def_var (dataFile, "m_k", NC_INT, 1, &nres_d, &mk_y);

  // PsiNk
  int PsiNk_y;
  err += nc_def_var (dataFile, "PsiN_k", NC_DOUBLE, 1, &nres_d, &PsiNk_y);

  // tau_Ak
  int tau_Ak_y;
  err += nc_def_var (dataFile, "tau_A", NC_DOUBLE, 1, &nres_d, &tau_Ak_y);

  // tau_Hk
  int tau_Hk_y;
  err += nc_def_var (dataFile, "tau_H", NC_DOUBLE, 1, &nres_d, &tau_Hk_y);

  // tau_Rk
  int tau_Rk_y;
  err += nc_def_var (dataFile, "tau_R", NC_DOUBLE, 1, &nres_d, &tau_Rk_y);

  // tau_Mk
  int tau_Mk_y;
  err += nc_def_var (dataFile, "tau_M", NC_DOUBLE, 1, &nres_d, &tau_Mk_y);

  // tau_thk
  int tau_thk_y;
  err += nc_def_var (dataFile, "tau_th", NC_DOUBLE, 1, &nres_d, &tau_thk_y);

  // tau_cxk
  int tau_cxk_y;
  err += nc_def_var (dataFile, "tau_cx", NC_DOUBLE, 1, &nres_d, &tau_cxk_y);

  // w_linear
  int w_linear_y;
  err += nc_def_var (dataFile, "w_linear", NC_DOUBLE, 1, &nres_d, &w_linear_y);

  // w_nonlinear
  int w_nonlinear_y;
  err += nc_def_var (dataFile, "w_nonlinear", NC_DOUBLE, 1, &nres_d, &w_nonlinear_y);

  // w_EB
  int w_EB_y;
  err += nc_def_var (dataFile, "w_EB", NC_DOUBLE, 1, &nres_d, &w_EB_y);

  // w_tk
  int w_tk_y;
  err += nc_def_var (dataFile, "w_phi_I", NC_DOUBLE, 1, &nres_d, &w_tk_y);

  // w_pk
  int w_pk_y;
  err += nc_def_var (dataFile, "w_theta_I", NC_DOUBLE, 1, &nres_d, &w_pk_y);

  // w_EB_0
  int w_EB0_y;
  err += nc_def_var (dataFile, "w_EB0", NC_DOUBLE, 1, &nres_d, &w_EB0_y);

  // w_EB_1
  int w_EB1_y;
  err += nc_def_var (dataFile, "w_EB1", NC_DOUBLE, 1, &nres_d, &w_EB1_y);

  // w_EB_2
  int w_EB2_y;
  err += nc_def_var (dataFile, "w_EB2", NC_DOUBLE, 1, &nres_d, &w_EB2_y);

  // w_nc_I0
  int w_nc_I0_y;
  err += nc_def_var (dataFile, "w_nc_I0", NC_DOUBLE, 1, &nres_d, &w_nc_I0_y);

  // w_nc_I1
  int w_nc_I1_y;
  err += nc_def_var (dataFile, "w_nc_I1", NC_DOUBLE, 1, &nres_d, &w_nc_I1_y);

  // w_nc_I2
  int w_nc_I2_y;
  err += nc_def_var (dataFile, "w_nc_I2", NC_DOUBLE, 1, &nres_d, &w_nc_I2_y);

  // w_pncI_0
  int w_pnc_I0_y;
  err += nc_def_var (dataFile, "w_pnc_I0", NC_DOUBLE, 1, &nres_d, &w_pnc_I0_y);

  // w_pncI_1
  int w_pnc_I1_y;
  err += nc_def_var (dataFile, "w_pnc_I1", NC_DOUBLE, 1, &nres_d, &w_pnc_I1_y);

  // w_pncI_2
  int w_pnc_I2_y;
  err += nc_def_var (dataFile, "w_pnc_I2", NC_DOUBLE, 1, &nres_d, &w_pnc_I2_y);

  // rho_sk
  int rho_sk_y;
  err += nc_def_var (dataFile, "rho_s", NC_DOUBLE, 1, &nres_d, &rho_sk_y);

  // delk
  int delk_y;
  err += nc_def_var (dataFile, "delta", NC_DOUBLE, 1, &nres_d, &delk_y);

  // rhothek
  int rhothek_y;
  err += nc_def_var (dataFile, "rho_theta_e", NC_DOUBLE, 1, &nres_d, &rhothek_y);

  // rhothik
  int rhothik_y;
  err += nc_def_var (dataFile, "rho_theta_i", NC_DOUBLE, 1, &nres_d, &rhothik_y);

  // WcritTek
  int WcritTek_y;
  err += nc_def_var (dataFile, "W_crit_Te", NC_DOUBLE, 1, &nres_d, &WcritTek_y);

  // WcritTik
  int WcritTik_y;
  err += nc_def_var (dataFile, "W_crit_Ti", NC_DOUBLE, 1, &nres_d, &WcritTik_y);

  // Wcritnek
  int Wcritnek_y;
  err += nc_def_var (dataFile, "W_crit_ne", NC_DOUBLE, 1, &nres_d, &Wcritnek_y);
 
  // w_ast_ek
  int w_ast_ek_y;
  err += nc_def_var (dataFile, "w_ast_e", NC_DOUBLE, 1, &nres_d, &w_ast_ek_y);

  // w_ast_ik
  int w_ast_ik_y;
  err += nc_def_var (dataFile, "w_ast_i", NC_DOUBLE, 1, &nres_d, &w_ast_ik_y);

  // w_ast_Ik
  int w_ast_Ik_y;
  err += nc_def_var (dataFile, "w_ast_I", NC_DOUBLE, 1, &nres_d, &w_ast_Ik_y);
  
  // Sk
  int Sk_y;
  err += nc_def_var (dataFile, "S", NC_DOUBLE, 1, &nres_d, &Sk_y);

  // tauk
  int tauk_y;
  err += nc_def_var (dataFile, "tau", NC_DOUBLE, 1, &nres_d, &tauk_y);

  // PEk
  int PEk_y;
  err += nc_def_var (dataFile, "P_perp", NC_DOUBLE, 1, &nres_d, &PEk_y);

  // PMk
  int PMk_y;
  err += nc_def_var (dataFile, "P_phi", NC_DOUBLE, 1, &nres_d, &PMk_y);

  // Dk
  int Dk_y;
  err += nc_def_var (dataFile, "D", NC_DOUBLE, 1, &nres_d, &Dk_y);

  // QEk
  int QEk_y;
  err += nc_def_var (dataFile, "Q_E", NC_DOUBLE, 1, &nres_d, &QEk_y);
   
  // Qek
  int Qek_y;
  err += nc_def_var (dataFile, "Q_e", NC_DOUBLE, 1, &nres_d, &Qek_y);
   
  // Qik
  int Qik_y;
  err += nc_def_var (dataFile, "Q_i", NC_DOUBLE, 1, &nres_d, &Qik_y);
  
  // alpbek
  int alpbek_y;
  err += nc_def_var (dataFile, "alpha_b_e", NC_DOUBLE, 1, &nres_d, &alpbek_y);

  // alpbik
  int alpbik_y;
  err += nc_def_var (dataFile, "alpha_b_i", NC_DOUBLE, 1, &nres_d, &alpbik_y);

  // alpck
  int alpck_y;
  err += nc_def_var (dataFile, "alpha_c",   NC_DOUBLE, 1, &nres_d, &alpck_y);

  // alppk
  int alppk_y;
  err += nc_def_var (dataFile, "alpha_p",   NC_DOUBLE, 1, &nres_d, &alppk_y);

  // Factor1
  int fac1_y;
  err += nc_def_var (dataFile, "Factor_1",  NC_DOUBLE, 1, &nres_d, &fac1_y);

  // Factor2
  int fac2_y;
  err += nc_def_var (dataFile, "Factor_2",  NC_DOUBLE, 1, &nres_d, &fac2_y);

  // Factor3
  int fac3_y;
  err += nc_def_var (dataFile, "Factor_3",  NC_DOUBLE, 1, &nres_d, &fac3_y);

  // Factor4
  int fac4_y;
  err += nc_def_var (dataFile, "Factor_4",  NC_DOUBLE, 1, &nres_d, &fac4_y);

  // Factor5
  int fac5_y;
  err += nc_def_var (dataFile, "Factor_5",  NC_DOUBLE, 1, &nres_d, &fac5_y);

  // Factor6
  int fac6_y;
  err += nc_def_var (dataFile, "Factor_6",  NC_DOUBLE, 1, &nres_d, &fac6_y);

  // Factor7
  int fac7_y;
  err += nc_def_var (dataFile, "Factor_7",  NC_DOUBLE, 1, &nres_d, &fac7_y);

  // Factor8
  int fac8_y;
  err += nc_def_var (dataFile, "Factor_8",  NC_DOUBLE, 1, &nres_d, &fac8_y);

  // Factor9
  int fac9_y;
  err += nc_def_var (dataFile, "Factor_9",  NC_DOUBLE, 1, &nres_d, &fac9_y);

  // Factor10
  int fac10_y;
  err += nc_def_var (dataFile, "Factor_10", NC_DOUBLE, 1, &nres_d, &fac10_y);

  // Factor11
  int fac11_y;
  err += nc_def_var (dataFile, "Factor_11", NC_DOUBLE, 1, &nres_d, &fac11_y);

  // Factor12
  int fac12_y;
  err += nc_def_var (dataFile, "Factor_12", NC_DOUBLE, 1, &nres_d, &fac12_y);

  // w_nc_i
  int w_nc_i_y;
  err += nc_def_var (dataFile, "w_nc_i", NC_DOUBLE, 1, &nres_d, &w_nc_i_y);
 
  // w_nc_i
  int w_nc_e_y;
  err += nc_def_var (dataFile, "w_nc_e", NC_DOUBLE, 1, &nres_d, &w_nc_e_y);
  
  // w_nc_i
  int w_th_i_y;
  err += nc_def_var (dataFile, "w_th_i", NC_DOUBLE, 1, &nres_d, &w_th_i_y);

  // w_para
  int w_para_y;
  err += nc_def_var (dataFile, "w_para", NC_DOUBLE, 1, &nres_d, &w_para_y);

  // w_psi
  int island_d, w_y;
  err += nc_def_dim (dataFile, "N_island", NISLAND, &island_d);
  err += nc_def_var (dataFile, "w_psi", NC_DOUBLE, 1, &island_d, &w_y);
  
  // rhs
  int rhs_d[2], rhs1_y, rhs2_y, rhs3_y, rhs1a_y, rhs2a_y, rhs3a_y;
  rhs_d[0] = nres_d;
  rhs_d[1] = island_d;
  err += nc_def_var (dataFile, "rhs1",  NC_DOUBLE, 2, rhs_d, &rhs1_y);
  err += nc_def_var (dataFile, "rhs2",  NC_DOUBLE, 2, rhs_d, &rhs2_y);
  err += nc_def_var (dataFile, "rhs3",  NC_DOUBLE, 2, rhs_d, &rhs3_y);
  err += nc_def_var (dataFile, "rhs1a", NC_DOUBLE, 2, rhs_d, &rhs1a_y);
  err += nc_def_var (dataFile, "rhs2a", NC_DOUBLE, 2, rhs_d, &rhs2a_y);
  err += nc_def_var (dataFile, "rhs3a", NC_DOUBLE, 2, rhs_d, &rhs3a_y);

  err += nc_enddef (dataFile);

  if (err != 0)
    {
      printf ("NEOCLASSICAL::WriteStage3Netcdfc: Error defining variables in Outputs/Stage3.nc\n");
      exit (1);
    }

  // Write data
  err += nc_put_var_double (dataFile, PsiN,          psi_x);
  err += nc_put_var_double (dataFile, n_e_y,         n_e_x);
  err += nc_put_var_double (dataFile, T_e_y,         T_e_x);
  err += nc_put_var_double (dataFile, n_i_y,         n_i_x);
  err += nc_put_var_double (dataFile, T_i_y,         T_i_x);
  err += nc_put_var_double (dataFile, n_I_y,         n_I_x);
  err += nc_put_var_double (dataFile, T_I_y,         T_I_x);
  err += nc_put_var_double (dataFile, Z_eff_y,       Z_eff_x);
  err += nc_put_var_double (dataFile, n_n_y,         n_n_x);
  err += nc_put_var_double (dataFile, w_t_y,         w_t_x);
  err += nc_put_var_double (dataFile, w_p_y,         w_p_x);
  err += nc_put_var_double (dataFile, w_E_y,         w_E_x);
  err += nc_put_var_double (dataFile, chip_y,        chip_x);
  err += nc_put_var_double (dataFile, chie_y,        chie_x);
  err += nc_put_var_double (dataFile, chin_y,        chin_x);
  err += nc_put_var_double (dataFile, chii_y,        chii_x);
  err += nc_put_var_int    (dataFile, mk_y,          mk_x);
  err += nc_put_var_double (dataFile, PsiNk_y,       PsiNk_x);
  err += nc_put_var_double (dataFile, tau_Ak_y,      tau_Ak_x);
  err += nc_put_var_double (dataFile, tau_Hk_y,      tau_Hk_x);
  err += nc_put_var_double (dataFile, tau_Rk_y,      tau_Rk_x);
  err += nc_put_var_double (dataFile, tau_Mk_y,      tau_Mk_x);
  err += nc_put_var_double (dataFile, tau_thk_y,     tau_thk_x);
  err += nc_put_var_double (dataFile, tau_cxk_y,     tau_cxk_x);
  err += nc_put_var_double (dataFile, w_linear_y,    w_linear_x);
  err += nc_put_var_double (dataFile, w_nonlinear_y, w_nonlinear_x);
  err += nc_put_var_double (dataFile, w_tk_y,        w_tk_x);
  err += nc_put_var_double (dataFile, w_pk_y,        w_pk_x);
  err += nc_put_var_double (dataFile, w_EB_y,        w_EB_x);
  err += nc_put_var_double (dataFile, w_EB0_y,       w_EB0_x);
  err += nc_put_var_double (dataFile, w_EB1_y,       w_EB1_x);
  err += nc_put_var_double (dataFile, w_EB2_y,       w_EB2_x);
  err += nc_put_var_double (dataFile, w_nc_I0_y,     w_nc_I0_x);
  err += nc_put_var_double (dataFile, w_nc_I1_y,     w_nc_I1_x);
  err += nc_put_var_double (dataFile, w_nc_I2_y,     w_nc_I2_x);
  err += nc_put_var_double (dataFile, w_pnc_I0_y,    w_pnc_I0_x);
  err += nc_put_var_double (dataFile, w_pnc_I1_y,    w_pnc_I1_x);
  err += nc_put_var_double (dataFile, w_pnc_I2_y,    w_pnc_I2_x);
  err += nc_put_var_double (dataFile, rho_sk_y,      rho_sk_x);
  err += nc_put_var_double (dataFile, delk_y,        delk_x);
  err += nc_put_var_double (dataFile, rhothek_y,     rhothek_x);
  err += nc_put_var_double (dataFile, rhothik_y,     rhothik_x);
  err += nc_put_var_double (dataFile, WcritTek_y,    WcritTek_x);
  err += nc_put_var_double (dataFile, WcritTik_y,    WcritTik_x);
  err += nc_put_var_double (dataFile, Wcritnek_y,    Wcritnek_x);
  err += nc_put_var_double (dataFile, w_ast_ek_y,    w_ast_ek_x);
  err += nc_put_var_double (dataFile, w_ast_ik_y,    w_ast_ik_x);
  err += nc_put_var_double (dataFile, w_ast_Ik_y,    w_ast_Ik_x);
  err += nc_put_var_double (dataFile, Sk_y,          Sk_x);
  err += nc_put_var_double (dataFile, tauk_y,        tauk_x);
  err += nc_put_var_double (dataFile, PEk_y,         PEk_x);
  err += nc_put_var_double (dataFile, PMk_y,         PMk_x);
  err += nc_put_var_double (dataFile, Dk_y,          Dk_x);
  err += nc_put_var_double (dataFile, QEk_y,         QEk_x);
  err += nc_put_var_double (dataFile, Qek_y,         Qek_x);
  err += nc_put_var_double (dataFile, Qik_y,         Qik_x);
  err += nc_put_var_double (dataFile, alpbek_y,      alpbek_x);
  err += nc_put_var_double (dataFile, alpbik_y,      alpbik_x);
  err += nc_put_var_double (dataFile, alpck_y,       alpck_x);
  err += nc_put_var_double (dataFile, alppk_y,       alppk_x);
  err += nc_put_var_double (dataFile, fac1_y,        fac1_x);
  err += nc_put_var_double (dataFile, fac2_y,        fac2_x);
  err += nc_put_var_double (dataFile, fac3_y,        fac3_x);
  err += nc_put_var_double (dataFile, fac4_y,        fac4_x);
  err += nc_put_var_double (dataFile, fac5_y,        fac5_x);
  err += nc_put_var_double (dataFile, fac6_y,        fac6_x);
  err += nc_put_var_double (dataFile, fac7_y,        fac7_x);
  err += nc_put_var_double (dataFile, fac8_y,        fac8_x);
  err += nc_put_var_double (dataFile, fac9_y,        fac9_x);
  err += nc_put_var_double (dataFile, fac10_y,       fac10_x);
  err += nc_put_var_double (dataFile, fac11_y,       fac11_x);
  err += nc_put_var_double (dataFile, fac12_y,       fac12_x);
  err += nc_put_var_double (dataFile, w_nc_i_y,      w_nc_i_x);
  err += nc_put_var_double (dataFile, w_nc_e_y,      w_nc_e_x);
  err += nc_put_var_double (dataFile, w_th_i_y,      w_th_i_x);
  err += nc_put_var_double (dataFile, w_para_y,      w_para_x);
  err += nc_put_var_double (dataFile, w_y,           w_x);
  err += nc_put_var_double (dataFile, rhs1_y,        rhs1_x);
  err += nc_put_var_double (dataFile, rhs2_y,        rhs2_x);
  err += nc_put_var_double (dataFile, rhs3_y,        rhs3_x);
  err += nc_put_var_double (dataFile, rhs1a_y,       rhs1a_x);
  err += nc_put_var_double (dataFile, rhs2a_y,       rhs2a_x);
  err += nc_put_var_double (dataFile, rhs3a_y,       rhs3a_x);
 
  if (err != 0)
    {
      printf ("NEOCLASSICAL::WriteStage3Netcdfc: Error writing Outputs/Stage3.nc\n");
      exit (1);
    }
  
  // Close file
  err += nc_close (dataFile);

  if (err != 0)
    {
      printf ("NEOCLASSICAL::WriteStage3Netcdfc: Error closing Outputs/Stage3.nc\n");
      exit (1);
    }

  // Clean up
  delete[] psi_x;      delete[] n_e_x;      delete[] T_e_x;      delete[] n_i_x;         delete[] T_i_x;        
  delete[] n_I_x;      delete[] T_I_x;      delete[] w_E_x;      delete[] Z_eff_x;       delete[] w_t_x;   
  delete[] w_p_x;      delete[] n_n_x;      delete[] chip_x;     delete[] chie_x;        delete[] chin_x;     
  delete[] chii_x;     delete[] tau_Hk_x;   delete[] tau_Rk_x;   delete[] tau_Mk_x;      delete[] tau_thk_x;  
  delete[] tau_cxk_x;  delete[] PsiNk_x;    delete[] w_linear_x; delete[] w_nonlinear_x; delete[] w_EB_x;     
  delete[] rho_sk_x;   delete[] delk_x;     delete[] rhothek_x;  delete[] rhothik_x;     delete[] WcritTek_x; 
  delete[] WcritTik_x; delete[] Wcritnek_x; delete[] w_ast_ek_x; delete[] w_ast_ik_x;    delete[] w_ast_Ik_x; 
  delete[] Sk_x;       delete[] alpbek_x;   delete[] alpbik_x;   delete[] alpck_x;       delete[] alppk_x;   
  delete[] fac1_x;     delete[] fac2_x;     delete[] fac3_x;     delete[] fac4_x;        delete[] fac5_x;     
  delete[] fac6_x;     delete[] fac7_x;     delete[] fac8_x;     delete[] fac9_x;        delete[] fac10_x;    
  delete[] fac11_x;    delete[] fac12_x;    delete[] w_nc_i_x;   delete[] w_nc_e_x;      delete[] w_th_i_x;  
  delete[] w_para_x;   delete[] w_EB0_x;    delete[] w_EB1_x;    delete[] w_EB2_x;       delete[] w_nc_I0_x;  
  delete[] w_nc_I1_x;  delete[] w_nc_I2_x;  delete[] w_pnc_I0_x; delete[] w_pnc_I1_x;    delete[] w_pnc_I2_x; 
  delete[] w_tk_x;     delete[] w_pk_x;     delete[] rhs1_x;     delete[] w_x;           delete[] mk_x;
  delete[] rhs2_x;     delete[] rhs3_x;     delete[] rhs1a_x;    delete[] rhs2a_x;       delete[] rhs3a_x;
  delete[] tauk_x;     delete[] PEk_x;      delete[] PMk_x;      delete[] Dk_x;          delete[] QEk_x;
  delete[] Qek_x;      delete[] Qik_x;      delete[] tau_Ak_x;
}

// #######################################################################################
// Calculate normalized quantities at rational surfaces for program PHASE and output nFile
// #######################################################################################
void Neoclassical::Get_Normalized ()
{
  FILE* file = OpenFilew ((char*) "Outputs/nFile");

  fprintf (file, "%4d %16.9e %16.9e\n", nres, tau_A, P0 /1.e19/e/1.e3);
 
  for (int j = 0; j < nres; j++)
    {
      double Sk  = tau_Rk (j) * Q_00 (j) /tau_A;
      double wkl = w_linear (j) * tau_A;
      double wke = w_EB (j) * tau_A;
      double wkn = w_nonlinear (j) * tau_A;
      double tm  = tau_Mk (j) /tau_A;
      double th  = tau_thk (j) /mu_00_i (j) /tau_A;
      double tx  = tau_cxk (j) /tau_A;
     
      printf ("m = %3d Psi = %10.3e r = %10.3e q = %10.3e rho = %10.3e a = %10.3e S = %10.3e tauM = %10.3e tauth = %10.3e taucx = %10.3e del = %10.3e\n",
	      mk (j), PsiNk (j), rk (j), qk (j), rhok (j), a /R_0, Sk, tm, th, tx, delk (j));

      fprintf (file, "%4d %4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	       mk (j),                      ntor,                        rk (j),                     qk (j),                     rhok (j),
	       a /R_0,                      Sk,                          tm,                         th,
	       sqrt (qk (j)/gk (j)/sk (j)), delk (j),                    wkl,                        wke,                        wkn, 
	       dnedrk (j) /1.e19,           dTedrk (j) /e/1.e3,          Wcritnek (j),               WcritTek (j),               WcritTik (j),
	       akk (j),                     gk (j),                      dPsidr (j),                 PsiNk (j),                  nek (j) /1.e19,
	       nik (j) /1.e19,              Tek (j) /e/1.e3,             Tik (j) /e/1.e3,            dnidrk (j) /1.e19,          dTidrk (j) /e/1.e3,
	       Factor1 (j) /1.e19/e/1.e3,   Factor2  (j) /1.e19/e/1.e3,  Factor3  (j) /1.e19/e/1.e3, Factor4  (j) /1.e19/e/1.e3,
	       Factor5 (j) /1.e19/e/1.e3,   Factor6  (j) /1.e19/e/1.e3,  Factor7  (j) /1.e19/e/1.e3, Factor8  (j) /1.e19/e/1.e3,
	       Factor9 (j) /1.e19/e/1.e3,   Factor10 (j) /1.e19/e/1.e3,  Factor11 (j) /1.e19/e/1.e3, Factor12 (j) /1.e19/e/1.e3, tx,
	       alpbek (j),                  alpbik (j),                  alpck (j),                  alppk (j),
	       rhothek (j),                 rhothik (j),                 eta_ek (j),                 eta_ik (j),                 chipk (j));
	       
    }
   fclose (file);
   
   printf ("\n");
}

// ######################
// Chandrasekhar function
// ######################
double Neoclassical::psi_fun (double x)
{
  return gsl_sf_erf (x) - (2./sqrt(M_PI)) * x * exp (-x*x);
}

 // ####################################
// Derivative of Chandrasekhar function
// ####################################
double Neoclassical::psi_fun_p (double x)
{
  return (2./sqrt(M_PI)) * x * exp (-x*x);
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Neoclassical::Rhs (double x, Array<double,1>& y, Array<double,1>& dydx)
{
  if (flag == 0)
    {
      double nuDi = (3.*sqrt(M_PI) /4.) * ((1. - 1. /2./x/x) * psi_fun (x) + psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /4.) * alphak (jj) * ((1. - x_iI (jj) * x_iI (jj) /2./x/x) * psi_fun (x /x_iI (jj)) + psi_fun_p (x /x_iI (jj))) /x;
      double nuEi = (3.*sqrt(M_PI) /2.) * (psi_fun (x) - psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /2.) * alphak (jj) * ((AI/AII) * psi_fun (x /x_iI (jj)) - psi_fun_p (x /x_iI (jj))) /x;
      double nuTi = 3. * nuDi + nuEi;
      
      double nuDI = (3.*sqrt(M_PI) /4.) * ((1. - 1. /2./x/x) * psi_fun (x) + psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /4.) * (1./alphak (jj)) * ((1. - x_Ii (jj) * x_Ii (jj) /2./x/x) * psi_fun (x /x_Ii (jj)) + psi_fun_p (x /x_Ii (jj))) /x;
      double nuEI = (3.*sqrt(M_PI) /2.) * (psi_fun (x) - psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /2.) * (1./alphak (jj)) * ((AII/AI) * psi_fun (x /x_Ii (jj)) - psi_fun_p (x /x_Ii (jj))) /x;
      double nuTI = 3. * nuDI + nuEI;
      
      double nuDe = (3.*sqrt(M_PI) /4.) * ((1. - 1. /2./x/x) * psi_fun (x) + psi_fun_p (x))
	+ (3.*sqrt(M_PI) /4.) * Zeffk (jj);
      double nuEe = (3.*sqrt(M_PI) /2.) * (psi_fun (x) - psi_fun_p (x));
      double nuTe = 3. * nuDe + nuEe;
      
      dydx (0) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 8.)  * nuDi / (x*x*x   + nu_P_i (jj) * nuDi) /(x*x*x   + nu_PS_i (jj) * nuTi);
      dydx (1) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 10.) * nuDi / (x*x*x   + nu_P_i (jj) * nuDi) /(x*x*x   + nu_PS_i (jj) * nuTi);
      dydx (2) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 12.) * nuDi / (x*x*x   + nu_P_i (jj) * nuDi) /(x*x*x   + nu_PS_i (jj) * nuTi);
      dydx (3) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 8.)  * nuDI / (x*x*x   + nu_P_I (jj) * nuDI) /(x*x*x   + nu_PS_I (jj) * nuTI);
      dydx (4) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 10.) * nuDI / (x*x*x   + nu_P_I (jj) * nuDI) /(x*x*x   + nu_PS_I (jj) * nuTI);
      dydx (5) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 12.) * nuDI / (x*x*x   + nu_P_I (jj) * nuDI) /(x*x*x   + nu_PS_I (jj) * nuTI);
      dydx (6) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 9.)  * nuDe / (x*x*x*x + nu_P_e (jj) * nuDe) /(x*x*x*x + nu_PS_e (jj) * nuTe);
      dydx (7) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 11.) * nuDe / (x*x*x*x + nu_P_e (jj) * nuDe) /(x*x*x*x + nu_PS_e (jj) * nuTe);
      dydx (8) = gt (jj) * (8./3./sqrt(M_PI)) * exp (-x*x) * pow (x, 13.) * nuDe / (x*x*x*x + nu_P_e (jj) * nuDe) /(x*x*x*x + nu_PS_e (jj) * nuTe);
    }
  else
    {
      double QE  = QEk  (jres);
      double Qe  = Qek  (jres);
      double Qi  = Qik  (jres);
      double tau = tauk (jres);
      double PM  = PMk  (jres);
      double PE  = PEk  (jres);
      double D   = Dk   (jres);
       
      gsl_complex Y1 = gsl_complex_rect (y (0), y (1));
      gsl_complex Y2 = gsl_complex_rect (y (4), y (5));
      
      gsl_complex I    = gsl_complex_rect (0., 1.);
      gsl_complex fac1 = gsl_complex_rect (- QE * (QE + Qi) + PM * PE * x*x*x*x,          (QE + Qi) * (PM + PE) * x*x);
      gsl_complex fac2 = gsl_complex_rect (PE * x*x + (1. + 1./tau) * PM * D*D * x*x*x*x, (QE + Qe) + (QE + Qi) * D*D * x*x);
      gsl_complex fac  = gsl_complex_div (fac1, fac2);
      gsl_complex rhs1 = gsl_complex_mul (fac, Y1);
      gsl_complex rhs2 = gsl_complex_mul (fac, Y2);
      rhs1             = gsl_complex_mul_real (rhs1, x*x);
      rhs2             = gsl_complex_mul_real (rhs2, x*x);

      dydx (0) = y (2);
      dydx (1) = y (3);
      dydx (2) = GSL_REAL (rhs1);
      dydx (3) = GSL_IMAG (rhs1);

      dydx (4) = y (6);
      dydx (5) = y (7);
      dydx (6) = GSL_REAL (rhs2);
      dydx (7) = GSL_IMAG (rhs2);
    }
}

// ########################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive fourth-order Runge-Kutta scheme
//
//   x       ... independent variable
//   y       ... array of dependent variables
//   h       ... step-length
//   t_err   ... actual truncation error per step 
//   acc     ... desired truncation error per step
//   S       ... step-length cannot change by more than this factor from
//                  step to step
//   rept    ... number of step recalculations		  
//   maxrept ... maximum allowable number of step recalculations		  
//   h_min   ... minimum allowable step-length
//   h_max   ... maximum allowable step-length
//   flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step whilst attempting to maintain 
//  constant truncation error per step of acc:
//
//   flag = 0 ... error is absolute
//   flag = 1 ... error is relative
//   flag = 2 ... error is mixed
//
//  If step-length falls below h_min then routine aborts
// #########################################################################
void Neoclassical::RK4Adaptive (double& x, Array<double,1>& y, double& h, double& t_err, double acc, double S, int& rept,
				int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
{
  int neqns = y.extent (0);
  Array<double,1> y0 (neqns), y1 (neqns);
  double hin = h;

  // Save initial data
  double x0 = x;
  y0        = y;

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
          err1  = fabs ((y(i) - y1(i)) /y(i));
          err2  = fabs (y(i) - y1(i));
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15) t_err = 1.e-15;

  // Calculate new step-length 
  double h_est = h * pow (fabs (acc /t_err), 0.2);

  // Prevent step-length from changing by more than factor S
  if (h_est /h > S)
    h *= S;
  else if (h_est /h < 1. /S)
    h /= S;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Abort if step-length falls below h_min
  if (fabs(h) < h_min)
    { 
      //printf ("Neoclassical::RK4Adpative: Warning - |h| < hmin at x = %9.2e\n", x);
      //exit (1);
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %9.2e hin = %9.2e err = %9.2e acc = %9.2e hest = %9.2e hout = %9.2e count = %3d\n", 
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
//  x  ... independent variable
//  y  ... array of dependent variables
//  h  ... step-length
// #####################################################################
void Neoclassical::RK4Fixed (double& x, Array<double,1>& y, double h)
{
  int neqns = y.extent(0);
  Array<double,1> dydx (neqns), k1 (neqns), k2 (neqns), k3 (neqns);
  Array<double,1> k4 (neqns), f (neqns);

  // Zeroth intermediate step 
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1 (i) = h * dydx (i);
      f  (i) = y (i) + k1 (i) /2.;
    }

  // First intermediate step 
  Rhs (x + h / 2., f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2 (i) = h * dydx (i);
      f  (i) = y (i) + k2 (i) /2.;
    }

  // Second intermediate step 
  Rhs (x + h / 2., f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3 (i) = h * dydx (i);
      f  (i) = y (i) + k3 (i);
    }

  // Third intermediate step 
  Rhs (x + h, f, dydx);
  for (int i = 0; i < neqns; i++)
    k4(i) = h * dydx (i);

  // Actual step 
  for (int i = 0; i < neqns; i++)
    y (i) += k1 (i) /6. + k2 (i) /3. + k3 (i) /3. + k4 (i) /6.;
  x += h;
}

// #################################
// Function to open file for writing
// #################################
FILE* Neoclassical::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("NEOCLASSICAL::OpenFilew: Error opening data-file %s\n", filename);
      exit (1);
    }
  return file;
}

// #############################
// Matrix multiplication routine
// #############################
void Neoclassical::Matrix_Mult (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i)
 {
   for (int is = 0; is < i; is++)
	{
	  for (int js = 0; js < i; js++)
	    {
	      double val = 0.;
	      
	      for (int ks = 0; ks < i; ks++)
		{
		  val += gsl_matrix_get (A, is, ks) * gsl_matrix_get (B, ks, js);
		}

	      gsl_matrix_set (C, is, js, val);
	    }
	}

 }

// #######################
// Matrix addition routine
// #######################
void Neoclassical::Matrix_Add (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i)
 {
   for (int is = 0; is < i; is++)
	{
	  for (int js = 0; js < i; js++)
	    {
	      double val = gsl_matrix_get (A, is, js) + gsl_matrix_get (B, is, js);

	      gsl_matrix_set (C, is, js, val);
	    }
	}

 }

// ##########################
// Matrix subtraction routine
// ##########################
void Neoclassical::Matrix_Sub (gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, int i)
 {
   for (int is = 0; is < i; is++)
	{
	  for (int js = 0; js < i; js++)
	    {
	      double val = gsl_matrix_get (A, is, js) - gsl_matrix_get (B, is, js);

	      gsl_matrix_set (C, is, js, val);
	    }
	}

 }

// #################################
// Function to open file for writing
// #################################
FILE* Neoclassical::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("NEOCLASSICAL::OpenFilew: Error opening data-file %s\n", filename);
      exit (1);
    }
  return file;
}

// ###################################
// Function to open file for appending
// ###################################
FILE* Neoclassical::OpenFilea (char* filename)
{
  FILE* file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("NEOCLASSICAL::OpenFilea: Error opening data-file %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to call operating system
// #################################
void Neoclassical::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("NEOCLASSICAL: Operating system call error executing %s\n", command);
      exit (1);
    }
}
