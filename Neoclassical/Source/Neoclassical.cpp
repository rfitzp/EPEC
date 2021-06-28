// Neoclassical.cpp

// PROGRAM ORGANIZATION:
// 
//        Neoclassical:: Neoclassical     ()
// void   Neoclassical:: Solve            (int _NEUTRAL, int _IMPURITY, int _EXB, int _INTP, int _INTF,
//			                   int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME, int _CATS)
// void   Neoclassical:: Read_Parameters  (int _NEUTRAL, int _IMPURITY, int _EXB, int _INTP, int _INTF,
// 				           int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME, int _CATS)
// void   Neoclassical:: Read_Equilibrium ()
// void   Neoclassical:: Read_Profiles    ()
// void   Neoclassical:: Get_Derived      ()
// void   Neoclassical:: Get_Viscosities  ()
// void   Neoclassical:: Get_Parameters   ()
// void   Neoclassical:: Get_Frequencies  ()
// void   Neoclassical:: Get_LayerWidths  ()
// void   Neoclassical:: Get_Normalized   ()
// double Neoclassical:: psi_fun          (double x)
// double Neoclassical:: psi_fun_p        (double x)
// void   Neoclassical:: RK4Adaptive      (double& x, Array<double,1>& y, double& h, 
//			                    double& t_err, double acc, double S, int& rept,
//			                    int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
// void   Neoclassical:: RK4Fixed         (double& x, Array<double,1>& y, double h)
// FILE*  Neoclassical:: OpenFilew        (char* filename)
// FILE*  Neoclassical:: OpenFilew        (char* filename)
// FILE*  Neoclassical:: OpenFilea        (char* filename)

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
void Neoclassical::Solve (int _NEUTRAL, int _IMPURITY, int _EXB, int _INTP, int _INTF,
			  int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME, int _CATS)
{
  // Read discharge parameters
  Read_Parameters (_NEUTRAL, _IMPURITY, _EXB, _INTP, _INTF, _INTC, _NTYPE, _NN, _LN, _YN, _TIME, _CATS);
  
  // Read equilibrium data
  Read_Equilibrium ();

  // Read profile data
  Read_Profiles ();

  // Calculate derived quantities at rational surface
  Get_Derived ();

  // Calculate neoclassical viscosities at rational surfaces
  Get_Viscosities ();

  // Calculate neoclassical parameters at rational surfaces
  Get_Parameters ();

  // Calculcate neoclassical frequencies at rational surfaces
  Get_Frequencies ();

  // Calculcate linear layer widths at rational surfaces
  Get_LayerWidths ();

  // Calculate normalized quantites at rational surface and output nFile
  Get_Normalized ();
}

// ####################################
// Read Neoclassical control parameters
// ####################################
void Neoclassical::Read_Parameters (int _NEUTRAL, int _IMPURITY, int _EXB, int _INTP, int _INTF,
				    int _INTC, int _NTYPE, double _NN, double _LN, double _YN, double _TIME, int _CATS)
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
 
  // Read namelist
  NameListRead (&IMPURITY, &NEUTRAL, &EXB, &INTP, &INTF, &INTC, &NTYPE, &NN, &LN, &SVN, &YN, &EN, &TIME, &COULOMB, &NSMOOTH, &CATS, &TAUMIN);

  // Override namelist values with command line option values
  if (_NEUTRAL > -1)
    NEUTRAL = _NEUTRAL;
  if (_IMPURITY > -1)
    IMPURITY = _IMPURITY;
  if (_EXB != -999999)
    EXB = _EXB;
  if (_TIME > 0.)
    TIME = _TIME;
  if (_NTYPE > -1)
    NTYPE = _NTYPE;
  if (_NN > 0.)
    NN = _NN;
  if (_LN > 0.)
    LN = _LN;
  if (_YN > 0.)
    YN = _YN;
  if (_INTP > -1)
    INTP = _INTP;
  if (_INTF > -1)
    INTF = _INTF;
  if (_INTC > -1)
    INTC = _INTC;
  if (_CATS > -1)
    CATS = _CATS;

  // Output input parameters
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Input parameters (from Inputs/Neoclassical.nml and command line options):\n");
  printf ("IMPURITY = %2d NEUTRAL = %2d EXB = %2d INTP = %2d INTF = %2d INTC = %2d NTYPE = %2d NN = %11.4e LN = %11.4e SVN = %11.4e YN = %11.4e EN = %11.4e TIME = %11.4e NSMOOTH = %3d CATS = %2d TAUMIN = %11.4e\n",
	  IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, NSMOOTH, CATS, TAUMIN);

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
  if (EXB < 0 || EXB > 1)
    {
      printf ("NEOCLASSICAL::Read_Parameters: Error invalid EXB value\n");
      exit (1);
    }
  if (TAUMIN <= -1.)
    {
      printf ("NEOCLASSICAL::Read_Parameters: Error invalid TAUMIN value\n");
      exit (1);
    }
   
  FILE* namelist = OpenFilew ((char*) "Inputs/InputParameters.txt");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "Input parameters (from Inputs/Neoclassical.nml and command line options):\n");
  fprintf (namelist, "IMPURITY = %2d NEUTRAL = %2d EXB = %2d INTP = %2d INTF = %2d INTC = %2d NTYPE = %2d NN = %11.4e LN = %11.4e SVN = %11.4e YN = %11.4e EN = %11.4e TIME = %11.4e NSMOOTH = %3d CATS = %2d TAUMIN = %11.4e\n",
	   IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, NSMOOTH, CATS, TAUMIN);
  fclose (namelist);
  
  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "Git Hash     = "); fprintf (monitor, GIT_HASH);     fprintf (monitor, "\n");
  fprintf (monitor, "Compile time = "); fprintf (monitor, COMPILE_TIME); fprintf (monitor, "\n");
  fprintf (monitor, "Git Branch   = "); fprintf (monitor, GIT_BRANCH);   fprintf (monitor, "\n\n");
  fprintf (monitor, "Input parameters (from Inputs/Neoclassical.nml and command line options):\n");
  fprintf (monitor, "IMPURITY = %2d NEUTRAL = %2d EXB = %2d INTP = %2d INTF = %2d INTC = %2d NTYPE = %2d NN = %11.4e LN = %11.4e SVN = %11.4e YN = %11.4e EN = %11.4e TIME = %11.4e NSMOOTH = %3d CATS = %2d TAUMIN = %11.4e\n",
	   IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, NSMOOTH, CATS, TAUMIN);
  fclose (monitor);
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
  else
    {
      CallSystem ("cd Inputs; rm -rf fFile; ln -s ../../Flux/Outputs/fFile fFile");
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

  for (int j = 0; j < nres; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mk (j), &rk (j), &sk (j), &gk (j), &gmk (j), &Ktk (j), &Kastk (j), &fck (j), &akk (j), &PsiNk (j),
		  &dPsidr (j), &Kthek (j), &in, &A2 (j), &q_hat (j), &C1 (j), &C2 (j)) != 17)
	{
	  printf ("NEOCLASSICAL:Read_Equilibrium: Error reading fFile (3)\n");
	  exit (1);
	}
      qk (j) = double (mk (j)) /double (ntor);
      printf ("m = %3d  r = %11.4e  s = %11.4e  g = %11.4e  gm = %11.4e  Kt = %11.4e  Kast = %11.4e  Kthe = %11.4e  fc = %11.4e  akk = %11.4e  PsiN = %11.4e  q_hat = %11.4e\n",
	      mk (j), rk (j), sk (j), gk (j), gmk (j), Ktk (j), Kastk (j), Kthek (j), fck (j), akk (j), PsiNk (j), q_hat (j));
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

  printf ("NI  = %3d  ZI  = %11.4e  AI  = %11.4e\n", NI,  ZI,  AI);
  printf ("NII = %3d  ZII = %11.4e  AII = %11.4e\n", NII, ZII, AII);

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
  alpha.resize  (NPSI);

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
  w_E    (0)  = wE.GetY    (0);
  w_t    (0)  = wt.GetY    (0);
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
  w_E    (NPSI-1) = wE.GetY    (wE.GetN()-1);
  w_t    (NPSI-1) = wt.GetY    (wt.GetN()-1);
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
      w_E    (j) = InterpolateField (wE, psi (j), 0);
      w_t    (j) = InterpolateField (wt, psi (j), 0);
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
  
  // Output profiles
  FILE* file = OpenFilew ((char*) "Outputs/Stage3/profiles.txt");
  for (int j = 0; j < NPSI; j++)
    if (psi (j) < PSIRAT)
      fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	       psi (j), rr (j),
	       n_e (j) /1.e19, dn_edr (j) /1.e19, T_e (j) /1.e3/e, dT_edr (j) /1.e3/e,
	       n_i (j) /1.e19, dn_idr (j) /1.e19, T_i (j) /1.e3/e, dT_idr (j) /1.e3/e,
	       n_I (j) /1.e19, dn_Idr (j) /1.e19, T_I (j) /1.e3/e, dT_Idr (j) /1.e3/e,
	       w_E (j) /1.e3, Quasi (j), Z_eff (j), alpha (j), n_b (j) /1.e19, dpsidr (j) /a, w_t (j) /1.e3, n_n (j) /1.e19,
	       dn_edP1 (j) /1.e19, dT_edP1 (j) /1.e3/e, dn_idP1 (j) /1.e19, dT_idP1 (j) /1.e3/e,
	       dn_edP2 (j) /1.e19, dT_edP2 (j) /1.e3/e, dn_idP2 (j) /1.e19, dT_idP2 (j) /1.e3/e,
	       dn_edP3 (j) /1.e19, dT_edP3 (j) /1.e3/e, dn_idP3 (j) /1.e19, dT_idP3 (j) /1.e3/e);
  fclose (file);
  
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

  // Output perpendicular diffusivities
  file = OpenFilew ((char*) "Outputs/Stage3/chip.txt");
  for (int j = 0; j < NPSI; j++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", psi (j), rr (j), chip (j), chie (j), chin (j), chii (j));
  fclose (file);
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
  chiik.resize  (nres);

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
      wEk    (j) = Interpolate (NPSI, rr, w_E,    rk (j), 0);
      wtk    (j) = Interpolate (NPSI, rr, w_t,    rk (j), 0);
      nbk    (j) = Interpolate (NPSI, rr, n_b,    rk (j), 0);
      Zeffk  (j) = Interpolate (NPSI, rr, Z_eff,  rk (j), 0);
      alphak (j) = Interpolate (NPSI, rr, alpha,  rk (j), 0);
      chipk  (j) = Interpolate (NPSI, rr, chip,   rk (j), 0);
      chiek  (j) = Interpolate (NPSI, rr, chie,   rk (j), 0);
      chink  (j) = Interpolate (NPSI, rr, chin,   rk (j), 0);
      chiik  (j) = Interpolate (NPSI, rr, chii,   rk (j), 0);
      rhok   (j) = (AI * (nik (j) + nbk (j)) + AII * nIk (j)) * m_p /rho0;

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
 
      printf ("m = %3d r = %10.3e ne = %10.3e Te = %10.3e ni = %10.3e Ti = %10.3e nI = %10.3e TI = %10.3e wE = %10.3e wt = %10.3e Z_eff = %10.3e NN = %10.3e rho = %10.3e chip = %10.3e chie = %10.3e chin = %10.3e chii = %10.3e\n",
	      mk(j), rk(j), nek(j)/1.e19, Tek(j)/e/1.e3, nik(j)/1.e19, Tik(j)/e/1.e3, nIk(j)/1.e19, TIk(j)/e/1.e3, wEk(j)/1.e3, wtk(j)/1.e3, Zeffk(j), NNk(j)/1.e19, rhok(j), chipk(j), chiek(j), chink(j), chii(j));
    }

  if (IMPURITY == 0)
    for (int j = 0; j < nres; j++)
      {
	Zeffk  (j) = ZI;
	alphak (j) = 1.e-15;
      }

  if (NEUTRAL == 0)
    for (int j = 0; j < nres; j++)
      NNk (j) = 0.;

  // Output factors
  FILE* filex = OpenFilew ((char*) "Outputs/Stage3/factor.txt");
  for (int j = 0; j < nres; j++)
    fprintf (filex, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     PsiNk (j),
	     Factor1 (j) /1.e19/e/1.e3, Factor2  (j) /1.e19/e/1.e3, Factor3  (j) /1.e19/e/1.e3, Factor4  (j) /1.e19/e/1.e3,
	     Factor5 (j) /1.e19/e/1.e3, Factor6  (j) /1.e19/e/1.e3, Factor7  (j) /1.e19/e/1.e3, Factor8  (j) /1.e19/e/1.e3,
	     Factor9 (j) /1.e19/e/1.e3, Factor10 (j) /1.e19/e/1.e3, Factor11 (j) /1.e19/e/1.e3, Factor12 (j) /1.e19/e/1.e3);
  fclose (filex);

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
      WcritTek (j) = pow (chiek (j) /v_T_ek (j) /(rk (j)*a) /(rk (j)*a/R_0) /sk (j) /double (ntor), 1./3.) * rk (j) * a;
      WcritTik (j) = pow (chiik (j) /v_T_ik (j) /(rk (j)*a) /(rk (j)*a/R_0) /sk (j) /double (ntor), 1./3.) * rk (j) * a;
      Wcritnek (j) = pow (chink (j) /v_T_ik (j) /(rk (j)*a) /(rk (j)*a/R_0) /sk (j) /double (ntor), 1./3.) * rk (j) * a;
    }

  // ------------------------------------------------------
  // Calculate diamagnetic frequencies at rational surfaces
  // ------------------------------------------------------
  eta_ek.resize   (nres); eta_ik.resize   (nres); eta_Ik.resize   (nres);
  w_ast_ek.resize (nres); w_ast_ik.resize (nres); w_ast_Ik.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      eta_ek   (j) = nek (j) * dTedrk (j) /dnedrk (j) /Tek (j);
      eta_ik   (j) = nik (j) * dTidrk (j) /dnidrk (j) /Tik (j);		      
      eta_Ik   (j) = nIk (j) * dTIdrk (j) /dnIdrk (j) /TIk (j);
      w_ast_ek (j) =   (qk (j) /gk (j)) * Tek (j) * dnedrk (j) * (1. + eta_ek (j)) /e      /nek (j) /(a*rk (j)) /fabs (B_0);
      w_ast_ik (j) = - (qk (j) /gk (j)) * Tik (j) * dnidrk (j) * (1. + eta_ik (j)) /e /ZI  /nik (j) /(a*rk (j)) /fabs (B_0);
      w_ast_Ik (j) = - (qk (j) /gk (j)) * TIk (j) * dnIdrk (j) * (1. + eta_Ik (j)) /e /ZII /nIk (j) /(a*rk (j)) /fabs (B_0);

      printf ("m = %3d  wE = %11.4e  w_ast_e = %11.4e  w_ast_i = %11.4e  w_ast_I = %11.4e\n",
	      mk (j), wEk (j) /1.e3, w_ast_ek (j) /1.e3, w_ast_ik (j) /1.e3, w_ast_Ik (j) /1.e3);
    }

  // Output diamagnetic frequencies
  FILE* file = OpenFilew ((char*) "Outputs/Stage3/diamagnetic.txt");
  for (int j = 0; j < nres; j++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     rk (j), wEk (j) /1.e3, w_ast_ek (j) /1.e3, w_ast_ik (j) /1.e3, w_ast_Ik(j) /1.e3, Zeffk(j), NNk(j) /1.e19, PsiNk(j));
  fclose (file);

  // -----------------------------------------
  // Calculate timescales at rational surfaces
  // -----------------------------------------
  rho_sk.resize (nres); tau_Hk.resize (nres); tau_Rk.resize (nres); tau_Mk.resize (nres); tau_thk.resize (nres); tau_cxk.resize (nres);
 

  for (int j = 0; j < nres; j++)
    {
      rho_sk  (j) = sqrt (AI * m_p * Tek (j)) /e /fabs (B_0) /gk (j);
      
      tau_Hk  (j) = tau_A * sqrt (rhok(j)) * R_0 /a /gk(j) /sk(j) /double (ntor);
      tau_Rk  (j) = mu_0 * a*a * rk(j)*rk(j) * nek (j) * e*e /nu_eek (j) /m_e;
      tau_Mk  (j) = a*a /chipk(j);
      tau_thk (j) = 1. /nu_iik (j) / (1. + (qk(j) * R_0 /rk(j) /a) * (qk(j) * R_0 /rk(j) /a) /akk(j));
      tau_cxk (j) = 1. /NNk(j) /SVN;

      printf ("m = %3d  P0 = %11.4e  tau_A = %11.4e  tau_R = %11.4e  tau_M = %11.4e  tau_th = %11.4e  tau_cx = %11.4e\n",
	      mk (j), P0 /1.e19/e/1.e3, tau_A, tau_Rk (j), tau_Mk (j), tau_thk (j), tau_cxk (j));
    }

  // Output timescales
  file = OpenFilew ((char*) "Outputs/Stage3/timescale.txt");
  for (int j = 0; j < nres; j++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     rk (j), log10 (tau_A), log10 (tau_Rk (j)), log10 (tau_Mk (j)), log10 (tau_thk (j)), log10 (tau_cxk (j)));
  fclose (file);

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

  // Output neoclassical collisionality data
  file = OpenFilew ((char*) "Outputs/Stage3/neoclassical.txt");
  for (int j = 0; j < nres; j++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", rk (j), gt (j),
	     log10 (nu_P_e (j)),  log10 (nu_P_i (j)),  log10 (nu_P_I (j)),
	     log10 (nu_PS_e (j)), log10 (nu_PS_i (j)), log10 (nu_PS_I (j)));
   fclose (file);
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

      printf ("m = %3d  mu_i = (%11.4e, %11.4e, %11.4e)  mu_I = (%11.4e, %11.4e, %11.4e)  mu_e = (%11.4e, %11.4e, %11.4e)\n",
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
  G_Ii_00.resize (nres); Q_00.resize    (nres); 
 
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

      for (int is = 0; is < 4; is++)
	{
	  for (int js = 0; js < 4; js++)
	    {
	      double val = 0.;
	      
	      for (int ks = 0; ks < 4; ks++)
		{
		  val += gsl_matrix_get (invH, is, ks) * gsl_matrix_get (G, ks, js);
		}

	      gsl_matrix_set (L, is, js, val);
	    }
	}

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
            
      gsl_matrix_free      (H);
      gsl_matrix_free      (G);
      gsl_matrix_free      (L);
      gsl_matrix_free      (invH);
      gsl_permutation_free (p);

      // .............................
      // Calculate electron parameters
      // .............................
      gsl_matrix* E = gsl_matrix_alloc (2, 2);
      gsl_matrix* Q = gsl_matrix_alloc (2, 2);

      gsl_matrix_set (E, 0, 0, Zeffk (j)                       + mu_00_e (j));
      gsl_matrix_set (E, 0, 1, (3./2.) * Zeffk (j)             + mu_01_e (j));
      gsl_matrix_set (E, 1, 0, (3./2.) * Zeffk (j)             + mu_01_e (j));
      gsl_matrix_set (E, 1, 1, sqrt(2.) + (13./4.) * Zeffk (j) + mu_11_e (j));

      gsl_permutation* pp = gsl_permutation_alloc (2);

      gsl_linalg_LU_decomp (E, pp, &s);

      gsl_linalg_LU_invert (E, pp, Q);

      Q_00 (j) = gsl_matrix_get (Q, 0, 0);

      gsl_matrix_free (E);
      gsl_matrix_free (Q);
      gsl_permutation_free (pp);
      
     printf ("m = %3d  L_ii = (%10.3e, %10.3e)  L_iI = (%10.3e, %10.3e)  L_Ii = (%10.3e, %10.3e)  L_II = (%10.3e, %10.3e)  G_ii = %10.3e  G_Ii = %10.3e  Q_00 = %10.3e\n",
	     mk (j), L_ii_00 (j),  L_ii_01 (j),  L_iI_00 (j),  L_iI_01 (j),  L_Ii_00 (j),  L_Ii_01 (j),  L_II_00 (j), L_II_01 (j), G_ii_00 (j), G_Ii_00 (j), Q_00 (j));
    }

  // Output parameters
  FILE* file = OpenFilew ((char*) "Outputs/Stage3/parameters.txt");
  for (int j = 0; j < nres; j++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     rk (j), eta_ik (j), L_ii_00 (j),  L_ii_01 (j),  L_iI_00 (j),  L_iI_01 (j),
	     L_Ii_00 (j),  L_Ii_01 (j),  L_II_00 (j), L_II_01 (j), G_ii_00 (j), G_Ii_00 (j), Q_00 (j));
  fclose (file);
}

// #####################################################################
// Calculate natural island propagation frequencies at rational surfaces
// #####################################################################
void Neoclassical::Get_Frequencies ()
{
  // .....................
  // Calculate frequencies
  // .....................
  w_linear.resize (nres); w_nonlinear.resize (nres); w_EB.resize   (nres);
  w_nc_ik.resize  (nres); w_nc_Ik.resize     (nres); w_E_Ik.resize (nres);

  for (int j = 0; j < nres; j++)
    {
      w_nc_Ik (j) = - L_II_00 (j) * w_ast_Ik (j)
		    - L_Ii_00 (j) * w_ast_ik (j)
		    + L_II_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j)
	            + L_Ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j);

      w_E_Ik (j) = (wtk (j) - w_ast_Ik (j) - Kthek (j) * w_nc_Ik (j)) /(1. - Kthek (j) * G_Ii_00 (j));

      if (EXB == 0)
	{
	  w_nc_ik (j) = - G_ii_00 (j) * wEk (j)
	                - L_ii_00 (j) * w_ast_ik (j)
	                - L_iI_00 (j) * w_ast_Ik (j)
	                + L_ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                + L_iI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  
	  w_linear    (j) = - double (ntor) * (wEk (j) + w_ast_ek (j));
	  w_nonlinear (j) = - double (ntor) * (wEk (j) + w_ast_ik (j) + w_nc_ik (j));
	  w_EB        (j) = - double (ntor) * (wEk (j));
	}
      else
	{
	  w_nc_ik (j) = - G_ii_00 (j) * w_E_Ik (j)
	                - L_ii_00 (j) * w_ast_ik (j)
	                - L_iI_00 (j) * w_ast_Ik (j)
	                + L_ii_01 (j) * (eta_ik (j) /(1. + eta_ik (j))) * w_ast_ik (j)
	                + L_iI_01 (j) * (eta_Ik (j) /(1. + eta_Ik (j))) * w_ast_Ik (j);
	  
	  w_linear    (j) = - double (ntor) * (w_E_Ik (j) + w_ast_ek (j));
	  w_nonlinear (j) = - double (ntor) * (w_E_Ik (j) + w_ast_ik (j) + w_nc_ik (j));
	  w_EB        (j) = - double (ntor) * (w_E_Ik (j));
	}

      printf ("m = %3d  w_linear = %11.4e  w_nonlinear = %11.4e  w_EB = %11.4e\n",
	      mk (j), w_linear (j) /1.e3, w_nonlinear (j) /1.e3, w_EB (j) /1.e3);
    }

  // ..................
  // Output frequencies
  // ..................
  FILE* file = OpenFilew ((char*) "Outputs/Stage3/omega.txt");
  for (int j = 0; j < nres; j++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     rk (j), w_linear (j) /1.e3, w_nonlinear (j) /1.e3, w_EB (j) /1.e3, PsiNk (j), wEk (j) /1.e3, w_E_Ik (j) /1.e3,
	     wtk (j) /1.e3, w_ast_Ik (j) /1.e3, Kthek (j) * (- G_Ii_00 (j) * w_E_Ik (j) + w_nc_Ik (j)) /1.e3);
  fclose (file);
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
      double tauE = a*a * rk (j)*rk (j) /(chink (j) + (2./3.) * chiek (j));
      double tauM = R_0*R_0 * qk (j)*qk (j) /chipk (j);

      Sk   (j) = Q_00 (j) * tau_Rk (j) /tau_Hk (j);
      tauk (j) = - w_ast_ik (j) /w_ast_ek (j);
      PEk  (j) = tau_Rk(j) /tauE;
      PMk  (j) = tau_Rk(j) /tauM;
      Dk   (j) = (5./3.) * pow (Sk (j), 1./3.) * rho_sk (j) /a/rk(j);
      QEk  (j) =   pow (Sk (j), 1./3.) * double (ntor) * wEk (j)      * tau_Hk (j);
      Qek  (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_ast_ek (j) * tau_Hk (j);
      Qik  (j) = - pow (Sk (j), 1./3.) * double (ntor) * w_ast_ik (j) * tau_Hk (j);

      // Prevent 1 + tauk from taking negative value
      if (tauk (j) < TAUMIN)
	tauk (j) = TAUMIN;

      printf ("m = %3d  tau = %11.4e  P_E = %11.4e  P_M = %11.4e  D = %11.4e  Q_E = %11.4e  Q_e = %11.4e  Q_i = %11.4e\n",
	    mk (j), tauk (j), PEk (j), PMk (j), Dk (j), QEk (j), Qek (j), Qik (j));
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
      double          C    = PEk (j) / (1. + tauk (j)) /Dk (j)/Dk (j);
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

// #######################################################################################
// Calculate normalized quantities at rational surfaces for program PHASE and output nFile
// #######################################################################################
void Neoclassical::Get_Normalized ()
{
  FILE* file = OpenFilew ((char*) "Outputs/nFile");

  fprintf (file, "%4d %16.9e %16.9e\n", nres, tau_A, P0 /1.e19/e/1.e3);
 
  for (int j = 0; j < nres; j++)
    {
      double dk  = M_PI * sqrt (double (ntor) * fabs (w_ast_ek (j))) * tau_Hk (j) * rk (j) * a
	/ (rho_sk (j) /rk(j) /a) /sqrt (tau_Rk (j) * Q_00 (j));
      double Sk  = tau_Rk (j) * Q_00 (j) /tau_A;
      double wkl = w_linear (j) * tau_A;
      double wke = w_EB (j) * tau_A;
      double wkn = w_nonlinear (j) * tau_A;
      double tm  = tau_Mk (j) /tau_A;
      double th  = tau_thk (j) /mu_00_i (j) /tau_A;
      double tx  = tau_cxk (j) /tau_A;

      printf ("m = %3d Psi = %10.3e r = %10.3e q = %10.3e rho = %10.3e a = %10.3e S = %10.3e tauM = %10.3e tauth = %10.3e taucx = %10.3e del_SCi = %10.3e del_true = %10.3e\n",
	      mk (j), PsiNk (j), rk (j), qk (j), rhok (j), a /R_0, Sk, tm, th, tx, dk, delk (j));

      fprintf (file, "%4d %4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	       mk (j),                      ntor,                        rk (j),                     qk (j),                    rhok (j),
	       a /R_0,                      Sk,                          tm,                         th,
	       sqrt (qk (j)/gk (j)/sk (j)), delk (j),                    wkl,                        wke,                       wkn, 
	       dnedrk (j) /1.e19,           dTedrk (j) /e/1.e3,          Wcritnek (j),               WcritTek (j),              WcritTik (j),
	       akk (j),                     gk (j),                      dPsidr (j),                 PsiNk (j),                 nek (j) /1.e19,
	       nik (j) /1.e19,              Tek (j) /e/1.e3,             Tik (j) /e/1.e3,            dnidrk (j) /1.e19,         dTidrk (j) /e/1.e3,
	       Factor1 (j) /1.e19/e/1.e3,   Factor2  (j) /1.e19/e/1.e3,  Factor3  (j) /1.e19/e/1.e3, Factor4  (j) /1.e19/e/1.e3,
	       Factor5 (j) /1.e19/e/1.e3,   Factor6  (j) /1.e19/e/1.e3,  Factor7  (j) /1.e19/e/1.e3, Factor8  (j) /1.e19/e/1.e3,
	       Factor9 (j) /1.e19/e/1.e3,   Factor10 (j) /1.e19/e/1.e3,  Factor11 (j) /1.e19/e/1.e3, Factor12 (j) /1.e19/e/1.e3, tx);
	       
      }
   fclose (file);

   if (INTP != 0)
     {
      char* filename = new char[MAXFILENAMELENGTH];
      sprintf (filename, "Outputs/nFiles/n.%d", int (TIME));
      file = OpenFilew (filename);
      fprintf (file, "%4d %16.9e %16.9e\n", nres, tau_A, P0 /1.e19/e/1.e3);
 
      for (int j = 0; j < nres; j++)
	{
	  double dk  = M_PI * sqrt (double (ntor) * fabs (w_ast_ek (j))) * tau_Hk (j) * rk (j) * a
	    / (rho_sk (j) /rk(j) /a) /sqrt (tau_Rk (j) * Q_00 (j));
	  double Sk  = tau_Rk (j) * Q_00 (j) /tau_A;
	  double wkl = w_linear (j) * tau_A;
	  double wke = w_EB (j) * tau_A;
	  double wkn = w_nonlinear (j) * tau_A;
	  double tm  = tau_Mk (j) /tau_A;
	  double th  = tau_thk (j) /mu_00_i (j) /tau_A;
	  double tx  = tau_cxk (j) /tau_A;

	  fprintf (file, "%4d %4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
		   mk (j),                      ntor,                        rk (j),                     qk (j),                     rhok (j),
		   a /R_0,                      Sk,                          tm,                         th,
		   sqrt (qk (j)/gk (j)/sk (j)), delk (j),                    wkl,                        wke,                        wkn, 
		   dnedrk (j) /1.e19,           dTedrk (j) /e/1.e3,          Wcritnek (j),               WcritTek (j),               WcritTik (j),
		   akk(j),                      gk (j),                      dPsidr (j),                 PsiNk (j),                  nek (j) /1.e19,
		   nik (j) /1.e19,              Tek (j) /e/1.e3,             Tik (j) /e/1.e3,            dnidrk (j) /1.e19,          dTidrk (j) /e/1.e3,
		   Factor1 (j) /1.e19/e/1.e3,   Factor2  (j) /1.e19/e/1.e3,  Factor3  (j) /1.e19/e/1.e3, Factor4  (j) /1.e19/e/1.e3,
		   Factor5 (j) /1.e19/e/1.e3,   Factor6  (j) /1.e19/e/1.e3,  Factor7  (j) /1.e19/e/1.e3, Factor8  (j) /1.e19/e/1.e3,
		   Factor9 (j) /1.e19/e/1.e3,   Factor10 (j) /1.e19/e/1.e3,  Factor11 (j) /1.e19/e/1.e3, Factor12 (j) /1.e19/e/1.e3 , tx);
	}
      fclose (file);

      sprintf (filename, "n.%d", int (TIME));
      file = OpenFilea ((char*) "Outputs/nFiles/Index");
      fprintf (file, "%s %19.6e\n", filename, TIME);
      fclose (file);

      file = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
      fprintf (file, "Wrote nFile %s\n", filename);
      fclose (file);
      
      delete[] filename;
     }
   
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
      double nuDi = (3.*sqrt(M_PI) /4.) * ((1. - 1. /2./x) * psi_fun (x) + psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /4.) * alphak (jj) * ((1. - x_iI (jj) /2./x) * psi_fun (x /x_iI (jj)) + psi_fun_p (x /x_iI (jj))) /x;
      double nuEi = (3.*sqrt(M_PI) /2.) * (psi_fun (x) - psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /2.) * alphak (jj) * ((AI/AII) * psi_fun (x /x_iI (jj)) - psi_fun_p (x /x_iI (jj))) /x;
      double nuTi = 3. * nuDi + nuEi;
      
      double nuDI = (3.*sqrt(M_PI) /4.) * ((1. - 1. /2./x) * psi_fun (x) + psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /4.) * (1./alphak (jj)) * ((1. - x_Ii (jj) /2./x) * psi_fun (x /x_Ii (jj)) + psi_fun_p (x /x_Ii (jj))) /x;
      double nuEI = (3.*sqrt(M_PI) /2.) * (psi_fun (x) - psi_fun_p (x)) /x
	+ (3.*sqrt(M_PI) /2.) * (1./alphak (jj)) * ((AII/AI) * psi_fun (x /x_Ii (jj)) - psi_fun_p (x /x_Ii (jj))) /x;
      double nuTI = 3. * nuDI + nuEI;
      
      double nuDe = (3.*sqrt(M_PI) /4.) * ((1. - 1. /2./x) * psi_fun (x) + psi_fun_p (x))
	+ (3.*sqrt(M_PI) /4.) * Zeffk (jj);
      double nuEe = (3.*sqrt(M_PI) /2.) * (psi_fun (x) - psi_fun_p (x));
      double nuTe = 3. * nuDe + nuEe;
      
      dydx (0) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 2.) * nuDi / (x   + nu_P_i (jj) * nuDi) /(x   + nu_PS_i (jj) * nuTi);
      dydx (1) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 3.) * nuDi / (x   + nu_P_i (jj) * nuDi) /(x   + nu_PS_i (jj) * nuTi);
      dydx (2) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 4.) * nuDi / (x   + nu_P_i (jj) * nuDi) /(x   + nu_PS_i (jj) * nuTi);
      dydx (3) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 2.) * nuDI / (x   + nu_P_I (jj) * nuDI) /(x   + nu_PS_I (jj) * nuTI);
      dydx (4) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 3.) * nuDI / (x   + nu_P_I (jj) * nuDI) /(x   + nu_PS_I (jj) * nuTI);
      dydx (5) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 4.) * nuDI / (x   + nu_P_I (jj) * nuDI) /(x   + nu_PS_I (jj) * nuTI);
      dydx (6) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 4.) * nuDe / (x*x + nu_P_e (jj) * nuDe) /(x*x + nu_PS_e (jj) * nuTe);
      dydx (7) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 5.) * nuDe / (x*x + nu_P_e (jj) * nuDe) /(x*x + nu_PS_e (jj) * nuTe);
      dydx (8) = gt (jj) * (4./3./sqrt(M_PI)) * exp (-x) * pow (x, 6.) * nuDe / (x*x + nu_P_e (jj) * nuDe) /(x*x + nu_PS_e (jj) * nuTe);
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
      gsl_complex fac1 = gsl_complex_rect (- QE * (QE - Qi) + PM * PE * x*x*x*x,       (QE - Qi) * (PM + PE) * x*x);
      gsl_complex fac2 = gsl_complex_rect (PE * x*x + (1. + tau) * PM * D*D * x*x*x*x, (QE - Qe) + (QE - Qi) * D*D * x*x);
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
// #########################################################################
void Neoclassical::RK4Adaptive (double& x, Array<double,1>& y, double& h, 
				double& t_err, double acc, double S, int& rept,
				int maxrept, double h_min, double h_max, int flag, 
				int diag, FILE* file)
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
      //printf ("Neoclassical::RK4Adpative: Warning - |h| < hmin at x = %11.4e\n", x);
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
