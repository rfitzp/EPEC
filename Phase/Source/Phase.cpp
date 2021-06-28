// Phase.cpp

// #####################
// PROGRAM ORGANIZATION:
// $####################

//        Phase:: Phase               ()
// void   Phase:: Solve               (int _STAGE5, int _INTF, int _INTN, int _INTU, int _OLD, int _LIN, int _MID, int _COPT,
// 		                       double _TSTART, double _TEND, double _SCALE, double _CHIR, double _IRMP, int _HIGH, int _RATS,
//                                     double _CORE, int _FREQ, double _FFAC)
// void   Phase:: Read_Data           (int _STAGE5, int _INTF, int _INTN, int _INTU, int _OLD, int _LIN, int _MID, int _COPT,
//                                     double _TSTART, double _TEND, double _SCALE, double _CHIR, double _IRMP, int _HIGH, int _RATS,
//                                     double _CORE, int _FREQ, double _FFAC)
// void   Phase:: Scan_Shift          ()
// void   Phase:: Calc_Velocity       ()
// void   Phase:: Initialize          ()
// void   Phase:: Save                ()
// void   Phase:: IslandDynamics      ()
// void   Phase:: CalcRMP             (double t)
// void   Phase:: CalcCoil            (double t, double& IU, double& IM, double& IL, double& PU, double& PM, double& PL)
// int    Phase:: Findk               ()
// double Phase:: FindMax             (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML)
// double Phase:: FindMin             (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML)
// void   Phase:: ThreeCoil           (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML, double Delta, double& fun, double& deriv, double& dderiv)
// void   Phase:: CalcChiZeta         (double t)
// void   Phase:: Pack                (Array<double,1> y)
// void   Phase:: Unpack              (Array<double,1> y)
// void   Phase:: PackRhs             (Array<double,1> XkRHS, Array<double,1> YkRHS,
//		                       Array<double,2> alphakpRHS, Array<double,2> betakpRHS, Array<double,1> dydt)
// double Phase:: GetNaturalFrequency (int j)
// double Phase:: GetActualFrequency  (int j)
// double Phase:: GetDeltaOmega       (int j)
// double Phase:: GetDeltaOmegaTheta  (int j)
// double Phase:: GetDeltaOmegaPhi    (int j)
// double Phase:: GetDeltaVPhi        (int j)
// double Phase:: GetDeltaVParallel   (int j)
// double Phase:: GetDeltaVEB         (int j)
// double Phase:: GetDeltaEr          (int j)
// void   Phase:: Rhs                 (double t, Array<double,1>& y, Array<double,1>& dydt)
// void   Phase:: RK4Adaptive         (double& x, Array<double,1>& y, double& h, 
//			               double& t_err, double acc, double S, int& rept,
//			               int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
// void   Phase:: RK4Fixed            (double& x, Array<double,1>& y, double h)
// FILE*  Phase:: OpenFilew           (char* filename)
// FILE*  Phase:: OpenFiler           (char* filename)
// FILE*  Phase:: OpenFilea           (char* filename)
// void   Phase:: CallSystem          (char* command)

#include "Phase.h"

// ###########
// Constructor
// ###########
Phase::Phase ()
{
  // -----------------------------------------------------
  // Set default values of adaptive integration parameters
  // -----------------------------------------------------
  acc      = 1.e-11;
  h0       = 1.e-2;
  hmin     = 1.e-2;
  hmax     = 1.e2;
  maxrept  = 50;
  omegamax = 1000.;

  // ----------------------
  // Set physical constants
  // ----------------------
  e         = GSL_CONST_MKSA_ELECTRON_CHARGE;
  epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
  mu_0      = GSL_CONST_MKSA_VACUUM_PERMEABILITY;
  m_p       = GSL_CONST_MKSA_MASS_PROTON;
  m_e       = GSL_CONST_MKSA_MASS_ELECTRON;
}

// ##############################
// Function to perform simulation
// ##############################
void Phase::Solve (int _STAGE5, int _INTF, int _INTN, int _INTU, int _NATS, int _OLD, int _LIN, int _MID, int _COPT,
		   double _TSTART, double _TEND, double _SCALE, double _CHIR, double _IRMP, int _HIGH, int _RATS,
		   double _CORE, int _FREQ, double _FFAC)
{
  // Read input data
  Read_Data (_STAGE5, _INTF, _INTN, _INTU, _NATS, _OLD, _LIN, _MID, _COPT, _TSTART, _TEND, _SCALE, _CHIR, _IRMP, _HIGH, _RATS,
	     _CORE, _FREQ, _FFAC);

  // Scan RMP phase shift
  Scan_Shift ();

  if (STAGE5 != 0)
    {
      // Calculate velocity factors
      Calc_Velocity ();
      
      // Perform island dynamics simulation
      IslandDynamics ();
    }

  // Clean up
  delete[] TCTRL; delete[] ICTRL; delete[] PCTRL;
  
  gsl_matrix_complex_free (EE);      gsl_matrix_complex_free (FF);
  gsl_vector_complex_free (EI);      gsl_vector_complex_free (EO);  
  gsl_vector_complex_free (DeltaU);  gsl_vector_complex_free (DeltaL);
  gsl_vector_complex_free (DeltaM);  gsl_vector_complex_free (ChiU);
  gsl_vector_complex_free (ChiL);    gsl_vector_complex_free (ChiM);
}

// ###########################
// Function to read input data
// ###########################
void Phase::Read_Data (int _STAGE5, int _INTF, int _INTN, int _INTU, int _NATS, int _OLD, int _LIN, int _MID, int _COPT,
		       double _TSTART, double _TEND, double _SCALE, double _CHIR, double _IRMP, int _HIGH, int _RATS,
		       double _CORE, int _FREQ, double _FFAC)
{
  // Output version information
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  
  // ......................................
  // Set default values of input parameters
  // ......................................
  NFLOW  = 200;

  STAGE5 = 1;

  INTF   = 0;
  INTN   = 0;
  INTU   = 0;
  NATS   = 0;
  OLD    = 0;

  LIN    = 0;
  FREQ   = 0;
  FFAC   = 0.;

  MID    = 2;
  IRMP   = -1.e6;
  IFLA   = 0;
  COPT   = 0;
  CORE   = 1.;

  DT     = 1.e-5;
  TSTART = 0.;
  TEND   = 1.e6;

  SCALE  = 2.;
  PMAX   = 4.;
  CHIR   = 1.;
  HIGH   = 1;
  RATS   = 1;
  
  // Read input data from namelists (Inputs/Phase.nml, Inputs/Waveform.nml)
  printf ("........................................................................................\n");
  printf ("Input parameters (from Inputs/Phase.nml, Inputs/Waveform.nml, and command line options):\n");
  printf ("........................................................................................\n");

  TCTRL = new double[MAXCONTROLPOINTNUMBER];
  ICTRL = new double[MAXCONTROLPOINTNUMBER];
  PCTRL = new double[MAXCONTROLPOINTNUMBER];

  NameListRead (&NFLOW, &STAGE5, &INTF, &INTN, &INTU, &NATS, &OLD, &FREQ, &LIN, &MID, &COPT,
		&DT, &TSTART, &TEND, &SCALE, &PMAX, &CHIR, &HIGH, &RATS, &CORE, &FFAC, 
		&NCTRL, TCTRL, ICTRL, PCTRL);

  TT.resize (NCTRL);
  
  // Override namelist values with command line options
  if (_STAGE5 > -1)
    STAGE5 = _STAGE5;
  if (_INTF > -1)
    INTF = _INTF;
  if (_INTN > -1)
    INTN = _INTN;
  if (_INTU > -1)
    INTU = _INTU;
  if (_OLD > -1)
    OLD = _OLD;
  if (_NATS > -1)
    NATS = _NATS;
  if (_LIN > -1)
    LIN = _LIN;
  if (_MID > -1)
    MID = _MID;
  if (_FREQ > -1)
    FREQ = _FREQ;
  if (_FFAC > -1.e6)
    FFAC = _FFAC;
  if (_HIGH > -1)
    HIGH = _HIGH;
  if (_RATS > -1)
    RATS = _RATS;
  if (_COPT > -1)
    COPT = _COPT;
  if (_CORE > -1.e6)
    CORE = _CORE;
  if (_TSTART > -1.e6)
    TSTART = _TSTART;
  if (_TEND > -1.e6)
    TEND = _TEND;
  if (_SCALE > -1.e6)
    SCALE = _SCALE;
  if (_CHIR > -1.e6)
    CHIR = _CHIR;
  if (_IRMP > -1.e6)
    {
      IFLA = 1;
      IRMP = _IRMP;
    }
  
  // .............................
  // Output calculation parameters
  // .............................
  printf ("NFLOW = %4d STAGE5 = %2d INTF = %2d INTN = %2d INTU = %2d OLD = %2d LIN = %2d MID = %2d COPT = %2d CORE = %11.4e HIGH = %2d RATS = %2d\n",
	  NFLOW, STAGE5, INTF, INTN, INTU, OLD, LIN, MID, COPT, CORE, HIGH, RATS);
  printf ("FREQ = %2d FFAC = %11.4e DT = %11.4e TSTART = %11.4e TEND = %11.4e SCALE = %11.4e PMAX = %11.4e CHIR = %11.4e NCTRL = %4d IRMP = %11.4e\n",
	  FREQ, FFAC, DT, TSTART, TEND, SCALE, PMAX, CHIR, NCTRL, IRMP);

  // ............
  // Sanity check
  // ............
  if (NFLOW <= 0)
    {
      printf ("PHASE:: NFLOW must be positive\n");
      exit (1);
    }
  if (DT < 0.)
    {
      printf ("PHASE:: DT must be positive\n");
      exit (1);
    }
  if (TEND < TSTART)
    {
      printf ("PHASE:: TEND must be greater than TSTART\n");
      exit (1);
    }
  if (CHIR < 0.9999999)
    {
      printf ("PHASE:: CHIR must be greater than unity\n");
      exit (1);
    }
  if (COPT < 0 || COPT > 3)
    {
      printf ("PHASE:: Invalid COPT value\n");
      exit (1);
    }
  if (FREQ < 0 || FREQ > 2)
    {
      printf ("PHASE:: Invalid FREQ value\n");
      exit (1);
    }
  if (MID < 1 || MID > 3)
    {
      printf ("PHASE:: Invalid MID value\n");
      exit (1);
    }

  for (int i = 0; i < NCTRL; i++)
    printf ("T = %11.4e  IRMP = %11.4e  PRMP/pi = %11.4e\n", TCTRL[i], ICTRL[i], PCTRL[i]/M_PI);

  FILE* namelist = OpenFilew ((char*) "Inputs/InputParameters.txt");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "Input parameters (from Inputs/Phase.nml and command line options):\n");
  fprintf (namelist, "NFLOW = %4d STAGE5 = %2d INTF = %2d INTN = %2d INTU = %2d NATS = %2d OLD = %2d LIN = %2d MID = %2d COPT = %2d CORE = %11.4e HIGH = %2d RATS = %2d \n",
	   NFLOW, STAGE5, INTF, INTN, INTU, NATS, OLD, LIN, MID, COPT, CORE, HIGH, RATS);
  fprintf (namelist, "FREQ = %2d FFAC = %11.4e DT = %11.4e TSTART = %11.4e TEND = %11.4e SCALE = %11.4e PMAX = %11.4e CHIR = %11.4e NCTRL = %4d IRMP = %11.4e\n",
	   FREQ, FFAC, DT, TSTART, TEND, SCALE, PMAX, CHIR, NCTRL, IRMP);
  fclose (namelist);
  
  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "Git Hash     = "); fprintf (monitor, GIT_HASH);     fprintf (monitor, "\n");
  fprintf (monitor, "Compile time = "); fprintf (monitor, COMPILE_TIME); fprintf (monitor, "\n");
  fprintf (monitor, "Git Branch   = "); fprintf (monitor, GIT_BRANCH);   fprintf (monitor, "\n\n");
  fprintf (monitor, "Input parameters (from Inputs/Phase.nml and command line options):\n");
  fprintf (monitor, "NFLOW = %4d STAGE5 = %2d INTF = %2d INTN = %2d INTU = %2d NATS = %2d OLD = %2d LIN = %2d MID = %2d COPT = %2d CORE = %11.4e HIGH = %2d RATS = %2d\n",
	   NFLOW, STAGE5, INTF, INTN, INTU, NATS, OLD, LIN, MID, COPT, CORE, HIGH, RATS);
  fprintf (monitor, "FREQ = %2d FFAC = %11.4e DT = %11.4e TSTART = %11.4e TEND = %11.4e SCALE = %11.4e PMAX = %11.4e CHIR = %11.4e NCTRL = %4d IRMP = %11.4e\n",
	   FREQ, FFAC, DT, TSTART, TEND, SCALE, PMAX, CHIR, NCTRL, IRMP);
  fclose (monitor);

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
      fFileInterp (fFileName, fFileTime, fFileNumber, TSTART);
    }
   else
    {
      CallSystem ("cd Inputs; rm -rf fFile; ln -s ../../Flux/Outputs/fFile fFile");
    }

  // ...........................
  // Read data from program FLUX
  // ...........................
  printf ("...............................\n");
  printf ("Reading data from program FLUX:\n");
  printf ("...............................\n");
  
  int ini; double inr; int NPSI;
  double Freal,  Fimag;
  double Ereal,  Eimag; 
  double EIreal, EIimag, EOreal, EOimag; 
 
  FILE* file = OpenFiler ((char*) "Inputs/fFile");
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &R_0, &B_0, &inr, &q95, &r95, &qrat, &rrat, &q0, &qa, &NPSI, &ini, &nres, &PSILIM, &PSIPED, &Pped, &PSIRAT) != 16)
    {
      printf ("PHASE::Error reading fFile (1)\n");
      exit (1);
    }

  for (int j = 0; j < NPSI; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &inr, &inr, &inr) != 3)
	{
	  printf ("NEOCLASSICAL: Error reading fFile (2)\n");
	  exit (1);
	}
    }

  A1.resize    (nres);
  qhatk.resize (nres);
  C1k.resize   (nres);
  C2k.resize   (nres);
  for (int j = 0; j < nres; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &ini, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &A1(j), &inr, &qhatk(j), &C1k(j), &C2k(j)) != 17)
	{
	  printf ("NEOCLASSICAL: Error reading fFile (3)\n");
	  exit (1);
	}
    }
  
  FF = gsl_matrix_complex_alloc (nres, nres);
  FFh.resize (nres, nres);

  for (int j = 0; j < nres; j++)
    for (int k = 0; k < nres; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal, &Fimag) != 4)
	 {
	   printf ("PHASE::Error reading fFile (4)\n");
	   exit (1);
	 }

	gsl_matrix_complex_set (FF, j, k, gsl_complex_rect (Freal, Fimag));
	FFh (j, k) = gsl_complex_abs (gsl_complex_rect (Freal, Fimag));
      }
  
  EE = gsl_matrix_complex_alloc (nres, nres);
  EEh.resize (nres, nres);
  xih.resize (nres, nres);
   
  for (int j = 0; j < nres; j++)
    for (int k = 0; k < nres; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal, &Eimag) != 4)
	 {
	   printf ("PHASE::Error reading fFile (5)\n");
	   exit (1);
	 }

	gsl_matrix_complex_set (EE, j, k, gsl_complex_rect (Ereal, Eimag));
	EEh (j, k) = gsl_complex_abs (gsl_complex_rect (Ereal, Eimag));
	xih (j, k) = - gsl_complex_arg (gsl_complex_rect (Ereal, Eimag));
      }
 
  EI = gsl_vector_complex_alloc (nres); EO = gsl_vector_complex_alloc (nres);
  
  for (int j = 0; j < nres; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf", &ini, &EIreal, &EIimag, &EOreal, &EOimag) != 5)
	{
	  printf ("PHASE::Error reading fFile (6)\n");
	  exit (1);
	}

      gsl_complex ei = gsl_complex_rect (EIreal, EIimag);
      gsl_complex eo = gsl_complex_rect (EOreal, EOimag);
      
      gsl_vector_complex_set (EI, j, ei);
      gsl_vector_complex_set (EO, j, eo);
    }
  
  fclose (file);

  epsi.resize (nres); sigi.resize (nres);
  epso.resize (nres); sigo.resize (nres);
  
  for (int j = 0; j < nres; j++)
    {
      epsi (j) =   gsl_complex_abs (gsl_vector_complex_get (EI, j));
      epso (j) =   gsl_complex_abs (gsl_vector_complex_get (EO, j));
      sigi (j) = - gsl_complex_arg (gsl_vector_complex_get (EI, j));
      sigo (j) = - gsl_complex_arg (gsl_vector_complex_get (EO, j));
    }

  printf ("R_0 = %11.4e  B_0 = %11.4e  nres = %3d\n", R_0, B_0, nres);

  printf ("E-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%9.2e,%9.2e) ", GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
      printf ("\n");
    }
  
  printf ("E-vectors:\n");
  for (int i = 0; i < nres; i++)
    printf ("EI = (%11.4e, %11.4e)  EO = (%11.4e, %11.4e)\n",
	    GSL_REAL (gsl_vector_complex_get (EI, i)), GSL_IMAG (gsl_vector_complex_get (EI, i)),
	    GSL_REAL (gsl_vector_complex_get (EO, i)), GSL_IMAG (gsl_vector_complex_get (EO, i)));

  // .................
  // Interpolate nFile
  // .................
  if (INTN != 0)
    {
      // Save pwd
      char pwd[MAXFILENAMELENGTH];
      getcwd (pwd, MAXFILENAMELENGTH);

      // Remove nFile
      CallSystem ("rm -rf Inputs/nFile");

      // Get nFiles directory
      char nFileDir[MAXFILENAMELENGTH];
      CallSystem ("greadlink -f Inputs/nFiles > nFileDir");
      FILE* nfd = OpenFiler ("nFileDir");
      fscanf (nfd, "%s", nFileDir);
      fclose (nfd);
      CallSystem ("rm nFileDir");
      
      // Read nFile data
      char           Basename[MAXFILENAMELENGTH];
      char           Filename[MAXFILENAMELENGTH];
      char           filename[MAXFILENAMELENGTH];
      vector<string> nFileName;
      vector<double> nFileTime;
      double         filetime;
      int            nFileNumber;
      
      printf ("Reading nFile data:\n");

      chdir (nFileDir);
      getcwd (Basename, MAXFILENAMELENGTH);
      strcat (Basename, "/");

      FILE* file = OpenFiler ((char*) "Index");

      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
          strcpy (Filename, Basename);
	  strcat (Filename, filename);

	  nFileName.push_back (Filename);
	  nFileTime.push_back (filetime);
	}
      nFileNumber = nFileName.size ();

      fclose (file);
      chdir (pwd);
 
      // Interpolate nFiles
      nFileInterp (nFileName, nFileTime, nFileNumber, TSTART);
    }
  else
    {
      CallSystem ("cd Inputs; rm -rf nFile; ln -s ../../Neoclassical/Outputs/nFile nFile");
    }
  
  // ...................................
  // Read data from program NEOCLASSICAL
  // ...................................
  printf (".......................................\n");
  printf ("Reading data from program NEOCLASSICAL:\n");
  printf (".......................................\n");
 
  file = OpenFiler ((char*) "Inputs/nFile");
  int nresn;
  if (fscanf (file, "%d %lf %lf", &nresn, &tau_A, &P0) != 3)
    {
      printf ("PHASE::Error reading nFile (1)\n");
      exit (1);
    }
  if (nresn != nres)
    {
      printf ("PHASE:: Warning - nresn != nres\n");
      if (nresn < nres)
	nres = nresn;
    }

  printf ("tau_A = %11.4e  P0 = %11.4e\n", tau_A, P0);

  // Normalize times
  for (int i = 0; i < NCTRL; i++)
    TT (i) = TCTRL[i] *1.e-3/tau_A;
  dTT    = DT         *1.e-3/tau_A;
  Tstart = TSTART     *1.e-3/tau_A;
  Tend   = TEND       *1.e-3/tau_A;

  mk.resize      (nres); ntor.resize     (nres); rk.resize       (nres); qk.resize       (nres); rhok.resize   (nres);
  a.resize       (nres); Sk.resize       (nres); taumk.resize    (nres); tautk.resize    (nres); tauxk.resize  (nres);
  fack.resize    (nres); delk.resize     (nres); wkl.resize      (nres); wke.resize      (nres); wkn.resize    (nres);
  dnedrk.resize  (nres); dTedrk.resize   (nres); Wcrnek.resize   (nres); WcrTek.resize   (nres); WcrTik.resize (nres);
  akk.resize     (nres); gk.resize       (nres); dPsiNdr.resize  (nres); PsiN.resize     (nres); nek.resize    (nres);
  nik.resize     (nres); Tek.resize      (nres); Tik.resize      (nres); dnidrk.resize   (nres); dTidrk.resize (nres);
  Factor1.resize (nres); Factor2.resize  (nres); Factor3.resize  (nres); Factor4.resize  (nres);
  Factor5.resize (nres); Factor6.resize  (nres); Factor7.resize  (nres); Factor8.resize  (nres);
  Factor9.resize (nres); Factor10.resize (nres); Factor11.resize (nres); Factor12.resize (nres);

  for (int j = 0; j < nres; j++)
    if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&mk      (j), &ntor     (j), &rk       (j), &qk       (j), &rhok   (j),
		&a       (j), &Sk       (j), &taumk    (j), &tautk    (j), &fack   (j),
		&delk    (j), &wkl      (j), &wke      (j), &wkn      (j),
		&dnedrk  (j), &dTedrk   (j), &Wcrnek   (j), &WcrTek   (j), &WcrTik (j),
		&akk     (j), &gk       (j), &dPsiNdr  (j), &PsiN     (j), &nek    (j),
		&nik     (j), &Tek      (j), &Tik      (j), &dnidrk   (j), &dTidrk (j),
		&Factor1 (j), &Factor2  (j), &Factor3  (j), &Factor4  (j),
		&Factor5 (j), &Factor6  (j), &Factor7  (j), &Factor8  (j),
		&Factor9 (j), &Factor10 (j), &Factor11 (j), &Factor12 (j), &tauxk (j)) != 42)
      {
	printf ("PHASE::Error reading nFile (2)\n");
	exit (1);
      }
  fclose (file);

  for (int j = 0; j < nres; j++)
    {
      if (taumk (j) < 0.)
	{
	  printf ("PHASE:: Error - taumk (%2d) < 0.\n", j);
	  exit (1);
	}
      if (tautk (j) < 0.)
	{
	  printf ("PHASE:: Error - tautk (%2d) < 0.\n", j);
	  exit (1);
	}
      if (Sk (j) < 0.)
	{
	  printf ("PHASE:: Error - Sk (%2d) < 0.\n", j);
	  exit (1);
	}
    }

  for (int j = 0; j < nres; j++)
    printf ("m = %3d h_r = %10.3e q = %10.3e g = %10.3e akk = %10.3e h_rho = %10.3e h_a = %10.3e S = %10.3e h_tauM = %10.3e h_tauth = %10.3e h_del = %10.3e A1 = %10.3e q_hat = %10.3e\n",
	    mk (j), rk (j), qk (j), gk (j), akk (j), rhok (j), a (j), Sk (j), taumk (j), tautk (j), delk (j) /(rk (j) * a (j) * R_0), A1 (j), qhatk (j));

  // Set Deltak+/- values
  Deltakp.resize (nres); Deltakm.resize (nres);

  Deltakm (0) = CHIR * PsiN (0);
  for (int j = 1; j < nres; j++)
    Deltakm (j) = CHIR * (PsiN (j) - PsiN (j-1));

  for (int j = 0; j < nres-1; j++)
    Deltakp (j) = CHIR * (PsiN (j+1) - PsiN (j));
  Deltakp (nres-1) = CHIR * (PSILIM - PsiN (nres-1));
  
  // ......................................
  // Interpolate uFiles, mFiles, and lFiles
  // ......................................
  if (INTU != 0)
    {
      // Save pwd
      char pwd[MAXFILENAMELENGTH];
      getcwd (pwd, MAXFILENAMELENGTH);

      if (MID == 3)
	{
	  // Remove mFile
	  CallSystem ("rm -rf Inputs/mFile");
	  
	  // Get mFiles directory
	  char mFileDir[MAXFILENAMELENGTH];
	  CallSystem ("greadlink -f Inputs/mFiles > mFileDir");
	  FILE* ufd = OpenFiler ("mFileDir");
	  fscanf (ufd, "%s", mFileDir);
	  fclose (ufd);
	  CallSystem ("rm mFileDir");
	  
	  // Read mFile data
	  char           Basename [MAXFILENAMELENGTH];
	  char           Filename [MAXFILENAMELENGTH];
	  char           mfilename[MAXFILENAMELENGTH];
	  vector<string> mFileName;
	  double         mfiletime;
	  vector<double> mFileTime;
	  int            mFileNumber;
	  
	  printf ("Reading mFile data:\n");
	  
	  chdir (mFileDir);
	  getcwd (Basename, MAXFILENAMELENGTH);
	  strcat (Basename, "/");
	  
	  file = OpenFiler ((char*) "Index");
	  
	  while (fscanf (file, "%s %lf", &mfilename, &mfiletime) == 2)
	    {
	      strcpy (Filename, Basename);
	      strcat (Filename, mfilename);
	      
	      mFileName.push_back (Filename);
	      mFileTime.push_back (mfiletime);
	    }
	  mFileNumber = mFileName.size ();
	  
	  fclose (file);
	  chdir (pwd);
	  
	  // Interpolate mFiles
	  mFileInterp (mFileName, mFileTime, mFileNumber, TSTART);
	}

      if (MID >= 2)
	{
	  // Remove uFile
	  CallSystem ("rm -rf Inputs/uFile");
	  
	  // Get uFiles directory
	  char uFileDir[MAXFILENAMELENGTH];
	  CallSystem ("greadlink -f Inputs/uFiles > uFileDir");
	  FILE* ufd = OpenFiler ("uFileDir");
	  fscanf (ufd, "%s", uFileDir);
	  fclose (ufd);
	  CallSystem ("rm uFileDir");
	  
	  // Read uFile data
	  char           Basename [MAXFILENAMELENGTH];
	  char           Filename [MAXFILENAMELENGTH];
	  char           ufilename[MAXFILENAMELENGTH];
	  vector<string> uFileName;
	  double         ufiletime;
	  vector<double> uFileTime;
	  int            uFileNumber;
	  
	  printf ("Reading uFile data:\n");
	  
	  chdir (uFileDir);
	  getcwd (Basename, MAXFILENAMELENGTH);
	  strcat (Basename, "/");
	  
	  file = OpenFiler ((char*) "Index");
	  
	  while (fscanf (file, "%s %lf", &ufilename, &ufiletime) == 2)
	    {
	      strcpy (Filename, Basename);
	      strcat (Filename, ufilename);
	      
	      uFileName.push_back (Filename);
	      uFileTime.push_back (ufiletime);
	    }
	  uFileNumber = uFileName.size ();
	  
	  fclose (file);
	  chdir (pwd);

	  // Interpolate uFiles
	  uFileInterp (uFileName, uFileTime, uFileNumber, TSTART);
	}

      if (MID >= 1)
	{
	  // Remove lFile
	  CallSystem ("rm -rf Inputs/lFile");
	  
	  // Get lFiles directory
	  char lFileDir[MAXFILENAMELENGTH];
	  CallSystem ("greadlink -f Inputs/lFiles > lFileDir");
	  FILE* lfd = OpenFiler ("lFileDir");
	  fscanf (lfd, "%s", lFileDir);
	  fclose (lfd);
	  CallSystem ("rm lFileDir");
	  
	  // Read lFile data
	  char           Basename [MAXFILENAMELENGTH];
	  char           Filename [MAXFILENAMELENGTH];
	  char           lfilename[MAXFILENAMELENGTH];
	  vector<string> lFileName;
	  double         lfiletime;
	  vector<double> lFileTime;
	  int            lFileNumber;
	  
	  printf ("Reading lFile data:\n");
	  
	  chdir (lFileDir);
	  getcwd (Basename, MAXFILENAMELENGTH);
	  strcat (Basename, "/");
	  
	  file = OpenFiler ((char*) "Index");
	  
	  while (fscanf (file, "%s %lf", &lfilename, &lfiletime) == 2)
	    {
	      strcpy (Filename, Basename);
	      strcat (Filename, lfilename);
	      
	      lFileName.push_back (Filename);
	      lFileTime.push_back (lfiletime);
	    }
	  lFileNumber = lFileName.size ();
	  
	  fclose (file);
	  chdir (pwd);
	  
	  // Interpolate lFiles
	  lFileInterp (lFileName, lFileTime, lFileNumber, TSTART);
	}
    }

  // ...........................
  // Read data from program GPEC
  // ...........................
  printf ("...............................\n");
  printf ("Reading data from program GPEC:\n");
  printf ("...............................\n");
  
  char    line[MAXULFILELINELENGTH]; char line1[MAXULFILELINELENGTH];
  double  v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;
  double* QIN = new double[nres];
  double* PSI = new double[nres];
  double* DRE = new double[nres];
  double* DIM = new double[nres];
  double* CRE = new double[nres];
  double* CIM = new double[nres];
  double* WWW = new double[nres];
  DeltaU = gsl_vector_complex_calloc (nres);
  DeltaL = gsl_vector_complex_calloc (nres);
  DeltaM = gsl_vector_complex_calloc (nres);
  ChiU   = gsl_vector_complex_calloc (nres);
  ChiL   = gsl_vector_complex_calloc (nres);
  ChiM   = gsl_vector_complex_calloc (nres);

  double SCALEFACTOR = SCALE;
  printf ("SCALEFACTOR = %11.4e\n", SCALEFACTOR);

  if (MID == 3)
    {
      printf ("Middle coil:\n");
      file = OpenFiler ((char*) "Inputs/mFile");
      
      for (int i = 0; i < 5; i++)
	fgets (line, MAXULFILELINELENGTH, file);
      fgets (line1, MAXULFILELINELENGTH, file);
      for (int i = 0; i < 2; i++)
	fgets (line, MAXULFILELINELENGTH, file);
      
      char* token = strtok (line1, " "); 
      token = strtok (NULL, " "); 
      token = strtok (NULL, " ");
      int nsingm = atoi (token);

      int nres1 = nres;
      if (nsingm < nres)
	{
	  printf ("PHASE:: Warning - nsingm < nres\n");
	  nres1 = nsingm;
	}
      if (nsingm > nres)
	{
	  printf ("PHASE:: Warning - nsingm > nres\n");
	  nres1 = nres;
	}
      
       for (int i = 0; i < nres1; i++)
	if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		    &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12) != 12)
	  {
	    printf ("Error reading mFile\n");
	    exit (1);
	  }
	else
	  {
	    QIN[i] = v1;
	    PSI[i] = v2;
	    DRE[i] = v9;
	    DIM[i] = v10;
	    WWW[i] = 2.*v11;
	    
	    DRE[i] /= 2.*M_PI * qk (i) /SCALEFACTOR/SCALEFACTOR;
	    DIM[i] /= 2.*M_PI * qk (i) /SCALEFACTOR/SCALEFACTOR;
	    
	    CRE[i] =   DIM[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    CIM[i] = - DRE[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    
	    gsl_vector_complex_set (DeltaM, i, gsl_complex_rect (DRE[i], DIM[i]));
	    gsl_vector_complex_set (ChiM,   i, gsl_complex_rect (CRE[i], CIM[i]));
	    
	    double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiM, i));
	    double WUNRE = 4. * sqrt (A1 (i) * Psi);
	    double WFULL = sqrt (FFh (i, i) * EEh (i, i)) * WUNRE;
	    
	    printf ("q = %11.4e  Psi = %11.4e  PsiN = %11.4e  Delta = (%11.4e, %11.4e)  Chi = (%11.4e, %11.4e)  W_UNRE = %11.4e  W_UNRE/W_GPEC = %11.4e  W_FULL/W_GPEC = %11.4e\n",
		    QIN[i], PSI[i], PsiN(i), DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
	  }
      fclose (file);

      if (fabs (QIN[0] - qk(0)) > 1.e-3)
	{
	  printf ("PHASE:: Error - minimum resonant q values do not match in nFile and mFile\n");
	  exit (1);
	}

      FILE* filem = OpenFilew ((char*) "Outputs/Stage4/mFile.txt");
      for (int i = 0; i < nres1; i++)
	{
	  double br_unrc = 1.e4 * EEh (i, i) * gsl_complex_abs (gsl_vector_complex_get (ChiM, i)) * (pow (rk (i), - double (mk (i))) - pow (rk (i), + double (mk (i)))) * fabs (B_0) /2./a (i);
	  double br_full = 1.e4 * EEh (i, i) * gsl_complex_abs (gsl_vector_complex_get (ChiM, i)) *  pow (rk (i), - double (mk (i)))                                    * fabs (B_0) /2./a (i);

	  fprintf (filem, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", QIN[i], PSI[i], WWW[i], br_unrc, br_full);
	}
      fclose (filem);
    }

  if (MID >= 2)
    {
      printf ("Upper coil:\n");
      file = OpenFiler ((char*) "Inputs/uFile");
      
      for (int i = 0; i < 5; i++)
	fgets (line, MAXULFILELINELENGTH, file);
      fgets (line1, MAXULFILELINELENGTH, file);
      for (int i = 0; i < 2; i++)
	fgets (line, MAXULFILELINELENGTH, file);
      
      char* token = strtok (line1, " "); 
      token = strtok (NULL, " "); 
      token = strtok (NULL, " ");
      int nsingu = atoi (token);
      
      int nres1 = nres;
      if (nsingu < nres)
	{
	  printf ("PHASE:: Warning - nsingu < nres\n");
	  nres1 = nsingu;
	}
      if (nsingu > nres)
	{
	  printf ("PHASE:: Warning - nsingu > nres\n");
	  nres1 = nres;
	}
      
      double sum = 0.;
      for (int i = 0; i < nres1; i++)
	if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		    &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12) != 12)
	  {
	    printf ("Error reading uFile\n");
	    exit (1);
	  }
	else
	  {
	    QIN[i] = v1;
	    PSI[i] = v2;
	    DRE[i] = v9;
	    DIM[i] = v10;
	    WWW[i] = 2.*v11;
	    
	    DRE[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
	    DIM[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
	    
	    CRE[i] =   DIM[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk (i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    CIM[i] = - DRE[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk (i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    
	    gsl_vector_complex_set (DeltaU, i, gsl_complex_rect (DRE[i], DIM[i]));
	    gsl_vector_complex_set (ChiU,   i, gsl_complex_rect (CRE[i], CIM[i]));
	    
	    double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiU, i));
	    double WUNRE = 4. * sqrt (A1 (i) * Psi);
	    double WFULL = sqrt (FFh (i, i) * EEh (i, i)) * WUNRE;
	    
	    printf ("q = %11.4e  Psi = %11.4e  PsiN = %11.4e  Delta = (%11.4e, %11.4e)  Chi = (%11.4e, %11.4e)  W_UNRE = %11.4e  W_UNRE/W_GPEC = %11.4e  W_FULL/W_GPEC = %11.4e\n",
		    QIN[i], PSI[i], PsiN(i), DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
	    
	  }
      fclose (file);

      if (fabs (QIN[0] - qk(0)) > 1.e-3)
	{
	  printf ("PHASE:: Error - minimum resonant q values do not match in nFile and uFile\n");
	  exit (1);
	}

      FILE* fileu = OpenFilew ((char*) "Outputs/Stage4/uFile.txt");
      for (int i = 0; i < nres1; i++)
	{
	  double br_unrc = 1.e4 * EEh (i, i) * gsl_complex_abs (gsl_vector_complex_get (ChiU, i)) * (pow (rk (i), - double (mk (i))) - pow (rk (i), + double (mk (i)))) * fabs (B_0) /2./a (i);
	  double br_full = 1.e4 * EEh (i, i) * gsl_complex_abs (gsl_vector_complex_get (ChiU, i)) *  pow (rk (i), - double (mk (i)))                                    * fabs (B_0) /2./a (i);

	  fprintf (fileu, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", QIN[i], PSI[i], WWW[i], br_unrc, br_full);
	}
      fclose (fileu);
    }

  if (MID >= 1)
    {
      printf ("Lower coil:\n");
      file = OpenFiler ((char*) "Inputs/lFile");
      
      for (int i = 0; i < 5; i++)
	fgets (line, MAXULFILELINELENGTH, file);
      fgets (line1, MAXULFILELINELENGTH, file);
      for (int i = 0; i < 2; i++)
	fgets (line, MAXULFILELINELENGTH, file);
      
      char* token = strtok (line1, " "); 
      token = strtok (NULL, " "); 
      token = strtok (NULL, " ");
      int nsingl = atoi (token);

      int nres1 = nres;
      if (nsingl < nres)
	{
	  printf ("PHASE:: Warning - nsingl < nres\n");
	  nres1 = nsingl;
	}
      if (nsingl > nres)
	{
	  printf ("PHASE:: Warning - nsingl > nres\n");
	  nres1 = nres;
	}
      
      for (int i = 0; i < nres1; i++)
	if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		    &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12) != 12)
	  {
	    printf ("Error reading lFile\n");
	    exit (1);
	  }
	else
	  {
	    QIN[i] = v1;
	    PSI[i] = v2;
	    DRE[i] = v9;
	    DIM[i] = v10;
	    WWW[i] = 2.*v11;
	    
	    DRE[i] /= 2.*M_PI * qk (i) /SCALEFACTOR/SCALEFACTOR;
	    DIM[i] /= 2.*M_PI * qk (i) /SCALEFACTOR/SCALEFACTOR;
	    
	    CRE[i] =   DIM[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    CIM[i] = - DRE[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    
	    gsl_vector_complex_set (DeltaL, i, gsl_complex_rect (DRE[i], DIM[i]));
	    gsl_vector_complex_set (ChiL,   i, gsl_complex_rect (CRE[i], CIM[i]));
	    
	    double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiL, i));
	    double WUNRE = 4. * sqrt (A1 (i) * Psi);
	    double WFULL = sqrt (FFh (i, i) * EEh (i, i)) * WUNRE;
	    
	    printf ("q = %11.4e  Psi = %11.4e  PsiN = %11.4e  Delta = (%11.4e, %11.4e)  Chi = (%11.4e, %11.4e)  W_UNRE = %11.4e  W_UNRE/W_GPEC = %11.4e  W_FULL/W_GPEC = %11.4e\n",
		    QIN[i], PSI[i], PsiN(i), DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
	  }
      fclose (file);

      if (fabs (QIN[0] - qk(0)) > 1.e-3)
	{
	  printf ("PHASE:: Error - minimum resonant q values do not match in nFile and lFile\n");
	  exit (1);
	}
      
      FILE* filel = OpenFilew ((char*) "Outputs/Stage4/lFile.txt");
      for (int i = 0; i < nres1; i++)
	{
	  double br_unrc = 1.e4 * EEh (i, i) * gsl_complex_abs (gsl_vector_complex_get (ChiL, i)) * (pow (rk (i), - double (mk (i))) - pow (rk (i), + double (mk (i)))) * fabs (B_0) /2./a (i);
	  double br_full = 1.e4 * EEh (i, i) * gsl_complex_abs (gsl_vector_complex_get (ChiL, i)) *  pow (rk (i), - double (mk (i)))                                    * fabs (B_0) /2./a (i);

	  fprintf (filel, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", QIN[i], PSI[i], WWW[i], br_unrc, br_full);
	}
      fclose (filel);
    } 
 
  file = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/q.txt");
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", q0, q95, qa, qrat, TSTART);
  fclose (file);

  file = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/Chi.txt");
  for (int i = 0; i < nres; i++)
    {
      if (MID == 3)
	fprintf (file, "T = %11.4e  m = %3d  ChiU = (%11.4e, %11.4e)  ChiM = (%11.4e, %11.4e)  ChiL = (%11.4e, %11.4e)\n",
		 TSTART, mk(i),
		 GSL_REAL (gsl_vector_complex_get (ChiU, i)), GSL_IMAG (gsl_vector_complex_get (ChiU, i)),
		 GSL_REAL (gsl_vector_complex_get (ChiM, i)), GSL_IMAG (gsl_vector_complex_get (ChiM, i)),
		 GSL_REAL (gsl_vector_complex_get (ChiL, i)), GSL_IMAG (gsl_vector_complex_get (ChiL, i)));
      else if (MID == 2)
	fprintf (file, "T = %11.4e  m = %3d  ChiU = (%11.4e, %11.4e)  ChiL = (%11.4e, %11.4e)\n",
		 TSTART, mk(i),
		 GSL_REAL (gsl_vector_complex_get (ChiU, i)), GSL_IMAG (gsl_vector_complex_get (ChiU, i)),
		 GSL_REAL (gsl_vector_complex_get (ChiL, i)), GSL_IMAG (gsl_vector_complex_get (ChiL, i)));
      else
	fprintf (file, "T = %11.4e  m = %3d  ChiL = (%11.4e, %11.4e)\n",
		 TSTART, mk(i),
		 GSL_REAL (gsl_vector_complex_get (ChiL, i)), GSL_IMAG (gsl_vector_complex_get (ChiL, i)));
    }
  fclose (file);

  delete[] QIN; delete[] DRE; delete[] DIM; delete[] CRE; delete[] CIM; delete[] WWW; delete[] PSI;
 }

// ########################################################################
// Function to calculate vacuum flux versus relative RMP coil current phase
// ########################################################################
void Phase::Scan_Shift ()
{
  FILE* file1 = OpenFilew ((char*) "Outputs/Stage4/chi.txt");
  FILE* file2 = OpenFilew ((char*) "Outputs/Stage4/zeta.txt");
  FILE* file3 = OpenFilew ((char*) "Outputs/Stage4/vac.txt");

  int    I   = 360;
  double one = 1.;
  for (int i = 0; i <= I; i++)
    {
      double      pha   = double (i) * PMAX*M_PI /double (I);
      gsl_complex eiku  = gsl_complex_polar (one, - pha);
      gsl_complex eikl  = gsl_complex_polar (one, + pha);
      gsl_complex eikuh = gsl_complex_polar (one, - pha/2.);
      gsl_complex eiklh = gsl_complex_polar (one, + pha/2.);

      fprintf (file1, "%e", pha /M_PI);
      fprintf (file2, "%e", pha /M_PI);
      fprintf (file3, "%e", pha /M_PI);
      for (int j = 0; j < nres; j++)
	{
	  gsl_complex hl, hu, hm, h;

	  if (MID == 3)
	    {
	      hl = gsl_vector_complex_get (ChiL, j);
	      hu = gsl_vector_complex_get (ChiU, j);
	      hm = gsl_vector_complex_get (ChiM, j);

	      hl = gsl_complex_mul (hl, eikl);
	      hu = gsl_complex_mul (hu, eiku);

	      h  = hl;
	      h  = gsl_complex_add (h, hu);
	      h  = gsl_complex_add (h, hm);
	    }
	  else if (MID == 2)
	    {
	      hl = gsl_vector_complex_get (ChiL, j);
	      hu = gsl_vector_complex_get (ChiU, j);

	      hl = gsl_complex_mul (hl, eiklh);
	      hu = gsl_complex_mul (hu, eikuh);

	      h  = hl;
	      h  = gsl_complex_add (h, hu);
	    }
	  else
	    {
	      h = gsl_vector_complex_get (ChiL, j);
	    }	  

	  double chi  =   gsl_complex_abs (h);
	  double zeta = - gsl_complex_arg (h);
	  double wv   = 4. * R_0 * fack (j) * sqrt (chi);
	    
	  fprintf (file1, " %e", chi);
	  fprintf (file2, " %e", zeta/M_PI);
	  fprintf (file3, " %e", wv);	  
	}
      fprintf (file1, "\n");
      fprintf (file2, "\n");
      fprintf (file3, "\n");
    }

  fclose (file1);
  fclose (file2);
  fclose (file3);

  FILE* file4 = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/vac.txt");

  double      pha   = PMAX*M_PI /2.;
  gsl_complex eiku  = gsl_complex_polar (one, - pha);
  gsl_complex eikl  = gsl_complex_polar (one, + pha);
  gsl_complex eikuh = gsl_complex_polar (one, - pha/2.);
  gsl_complex eiklh = gsl_complex_polar (one, + pha/2.);
  
  for (int j = 0; j < nres; j++)
    {
      gsl_complex hl, hu, hm, h;

      if (MID == 3)
	{
	  hl = gsl_vector_complex_get (ChiL, j);
	  hu = gsl_vector_complex_get (ChiU, j);
	  hm = gsl_vector_complex_get (ChiM, j);
	  
	  hl = gsl_complex_mul (hl, eikl);
	  hu = gsl_complex_mul (hu, eiku);
	  
	  h  = hl;
	  h  = gsl_complex_add (h, hu);
	  h  = gsl_complex_add (h, hm);
	}
      else if (MID == 2)
	{
	  hl = gsl_vector_complex_get (ChiL, j);
	  hu = gsl_vector_complex_get (ChiU, j);

	  hl = gsl_complex_mul (hl, eiklh);
	  hu = gsl_complex_mul (hu, eikuh);
	  
	  h  = hl;
	  h  = gsl_complex_add (h, hu);
	}
      else
	{
	  h = gsl_vector_complex_get (ChiL, j);
	}	  
      
      double chi     =   gsl_complex_abs (h);
      double zeta    = - gsl_complex_arg (h);
      double wv      = 4. * fack (j) * sqrt (chi) /a (j);
      double br_unrc = 1.e4 * EEh (j, j) * chi * (pow (rk (j), - double (mk (j))) - pow (rk (j), + double (mk (j)))) * fabs (B_0) /2./a (j);
      double br_full = 1.e4 * EEh (j, j) * chi *  pow (rk (j), - double (mk (j)))                                    * fabs (B_0) /2./a (j);
      
      fprintf (file4, "%3d %16.9e %16.9e %16.9e %16.9e %16.9e\n", mk (j), rk (j), wv, TSTART, br_unrc, br_full);
    }
  fclose (file4);
}

// ######################################
// Function to calculate velocity factors
// ######################################
void Phase::Calc_Velocity ()
{
  j0p.resize  (NFLOW);
  j1p.resize  (NFLOW);
  torp.resize (nres, NFLOW);
  tort.resize (nres, NFLOW);
  natp.resize (nres, nres, NFLOW);
  natt.resize (nres, nres, NFLOW);

  for (int i = 0; i < NFLOW; i++)
    {
      j0p (i) = gsl_sf_bessel_zero_J0 (i+1);
      j1p (i) = gsl_sf_bessel_zero_J1 (i+1);
    }

  for (int j = 0; j < nres; j++)
    for (int i = 0; i < NFLOW; i++)
      {
	double yp = gsl_sf_bessel_J1 (j1p (i) * rk (j)) /rk (j);
	double zp = gsl_sf_bessel_J0 (j0p (i) * rk (j));
	double J1 = gsl_sf_bessel_J1 (j0p (i));
	double J2 = gsl_sf_bessel_Jn (2, j1p (i));

	torp (j, i) = double (mk   (j) * mk   (j)) * yp*yp /rhok (j) /a (j)/a (j) /J2/J2;
	tort (j, i) = double (ntor (j) * ntor (j)) * zp*zp /rhok (j)              /J1/J1;
      }

  for (int j = 0; j < nres; j++)
    for (int k = 0; k < nres; k++)
      for (int i = 0; i  < NFLOW; i++)
	{
	  double ypj = gsl_sf_bessel_J1 (j1p (i) * rk (j)) /rk (j);
	  double zpj = gsl_sf_bessel_J0 (j0p (i) * rk (j));
	  double ypk = gsl_sf_bessel_J1 (j1p (i) * rk (k)) /rk (k);
	  double zpk = gsl_sf_bessel_J0 (j0p (i) * rk (k));
	  
	  natp (j, k, i) = double (mk (j)) * ypj /double (mk (k)) /ypk;
	  natt (j, k, i) = zpj /zpk;
	}			     
}

// ##################################################
// Function to initialize island dynamics calculation
// ##################################################
void Phase::Initialize ()
{
  // Initialize variables
  for (int j = 0; j < nres; j++)
    {
      Psik (j) = 0.;
      phik (j) = 0.;
      Xk   (j) = 0.;
      Yk   (j) = 0.;
 
      for (int i = 0; i < NFLOW; i++)
	{
	  alphakp (j, i) = 0.;
	  betakp  (j, i) = 0.;
	}
    }

  for (int j = 0; j < nres; j++)
    {
      lock (j) = 0;
      ww   (j) = GetActualFrequency (j);
    }

  // Reinitialize from previous run
  if (OLD == 1)
    {
      printf ("Loading previous calculation from file Outputs/sFile\n");
      
      FILE* file = OpenFiler ((char*) "Outputs/sFile");

      int _nres, _NFLOW;

      if (fscanf (file, "%d %d", &_nres, &_NFLOW) != 2)
	{
	  printf ("PHASE: Error reading sFile (1)\n");
	  exit (1);
	}
      
      Array<double,1> _Psik    (_nres);
      Array<double,1> _phik    (_nres);
      Array<double,2> _alphakp (_nres, _NFLOW);
      Array<double,2> _betakp  (_nres, _NFLOW);
      Array<int,1>    _lock    (_nres);
      Array<double,1> _ww      (_nres);

      int in;
      for (int j = 0; j < _nres; j++)
	if (fscanf (file, "%d %lf %lf %d %lf\n", &in, &_Psik (j), &_phik (j), &_lock (j), &_ww (j)) != 5)
	  {
	    printf ("PHASE: Error reading sFile (2)\n");
	    exit (1);
	  }

      for (int j = 0; j < _nres; j++)
	for (int i = 0; i < _NFLOW; i++)
	  if (fscanf (file, "%d %d %lf %lf\n", &in, &in, &_alphakp (j, i), &_betakp (j, i)) != 4)
	    {
	      printf ("PHASE: Error reading sFile (3)\n");
	      exit (1);
	    }
        
      fclose (file);

      if (_nres > nres)
	_nres = nres;
      if (_NFLOW >  NFLOW)
	_NFLOW = NFLOW;

      for (int j = 0; j < _nres; j++)
	{
	  Psik (j) = _Psik (j);
	  phik (j) = _phik (j);
	  Xk   (j) = Psik (j) * cos (phik (j));
	  Yk   (j) = Psik (j) * sin (phik (j));
 	  lock (j) = _lock (j);
	  ww   (j) = _ww   (j);
	}

      for (int j = _nres; j < nres; j++)
	{
	  Psik (j) = 0.;
	  phik (j) = 0.;
	  Xk   (j) = 0.;
	  Yk   (j) = 0.;
 	  lock (j) = 0.;
	  ww   (j) = GetActualFrequency (j);
	}

      for (int j = 0; j < _nres; j++)
	{
	  for (int i = 0; i < _NFLOW; i++)
	    {
	      alphakp (j, i) = _alphakp (j, i);
	      betakp  (j, i) = _betakp  (j, i);
	    }
	  for (int i = _NFLOW; i < NFLOW; i++)
	    {
	      alphakp (j, i) = 0.;
	      betakp  (j, i) = 0.;
	    }
	}

      for (int j = _nres; j < nres; j++)
	{
	  for (int i = 0; i < NFLOW; i++)
	    {
	      alphakp (j, i) = 0.;
	      betakp  (j, i) = 0.;
	    }
	}
    }
}

// ############################################
// Function to save island dynamics calculation
// ############################################
void Phase::Save ()
{
  printf ("Saving calculation in file Outputs/sFile\n");
  
  FILE* file = OpenFilew ("Outputs/sFile");

  fprintf (file, "%d %d\n", nres, NFLOW);

  for (int j = 0; j < nres; j++)
    fprintf (file, "%3d %19.6e %19.6e %2d %19.6e\n", j, Psik (j), phik (j), lock (j), ww (j));

  for (int j = 0; j < nres; j++)
    for (int i = 0; i < NFLOW; i++)
      fprintf (file, "%3d %4d %19.6e %19.6e\n", j, i, alphakp (j, i), betakp (j, i));
  
  fclose (file);
}

// ##############################################
// Function to perform island dynamics simulation
// ##############################################
void Phase::IslandDynamics ()
{
  // Initialize calculation
  Psik.resize    (nres);
  phik.resize    (nres);
  Xk.resize      (nres);
  Yk.resize      (nres);
  alphakp.resize (nres, NFLOW);
  betakp.resize  (nres, NFLOW);
  lock.resize    (nres);
  ww.resize      (nres);

  Initialize ();

  // Integrate equations of motion
  chi.resize  (nres);
  zeta.resize (nres);

  double          t, h, t_err;
  int             rept, step = 0, skip = 0; count = 0;
  Array<double,1> y    (2*nres*(1+NFLOW));
  Array<double,1> dydt (2*nres*(1+NFLOW));
  Array<double,1> wwo  (nres);

  t = Tstart;
  h = h0;
  Pack (y);
 
  int    cnt    = 0;
  double dt     = 0.; 
  FILE*  file   = OpenFilew ((char*) "Outputs/Stage5/ntor.txt");
  FILE*  file1  = OpenFilew ((char*) "Outputs/Stage5/Psi.txt");
  FILE*  file2  = OpenFilew ((char*) "Outputs/Stage5/phi.txt");
  FILE*  file3  = OpenFilew ((char*) "Outputs/Stage5/W.txt");
  FILE*  file3a = OpenFilew ((char*) "Outputs/Stage5/Wdel.txt");
  FILE*  file4  = OpenFilew ((char*) "Outputs/Stage5/omega.txt");
  FILE*  file4a = OpenFilew ((char*) "Outputs/Stage5/omega0.txt");
  FILE*  file4b = OpenFilew ((char*) "Outputs/Stage5/Omegat.txt");
  FILE*  file4c = OpenFilew ((char*) "Outputs/Stage5/Omegap.txt");
  FILE*  file4d = OpenFilew ((char*) "Outputs/Stage5/Er.txt");
  FILE*  file5  = OpenFilew ((char*) "Outputs/Stage5/rmp.txt");
  FILE*  file6  = OpenFilew ((char*) "Outputs/Stage5/mirnov.txt");
  FILE*  file8  = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/scan.txt");  
  FILE*  file9  = OpenFilew ((char*) "Outputs/Stage5/deltane.txt");
  FILE*  file10 = OpenFilew ((char*) "Outputs/Stage5/deltaTe.txt");
  FILE*  file11 = OpenFilew ((char*) "Outputs/Stage5/dne.txt");
  FILE*  file12 = OpenFilew ((char*) "Outputs/Stage5/dTe.txt");
  FILE*  file13 = OpenFilew ((char*) "Outputs/Stage5/rkminus.txt");
  FILE*  file14 = OpenFilew ((char*) "Outputs/Stage5/rkplus.txt");
  FILE*  file15 = OpenFilew ((char*) "Outputs/Stage5/rkminusne.txt");
  FILE*  file16 = OpenFilew ((char*) "Outputs/Stage5/rkplusne.txt");
  FILE*  file17 = OpenFilew ((char*) "Outputs/Stage5/rkminusTe.txt");
  FILE*  file18 = OpenFilew ((char*) "Outputs/Stage5/rkplusTe.txt");
  FILE*  file19 = OpenFilew ((char*) "Outputs/Stage5/results.txt");
  FILE*  file20 = OpenFilew ((char*) "Outputs/Stage5/optimize.txt");

  fprintf (file, "%d\n", ntor (0));
  fclose (file);
  
  printf ("......................\n");
  printf ("Performing simulation:\n");
  printf ("......................\n");
  do
    {
      for (int j = 0; j < nres; j++)
	wwo (j) = ww (j);
      
      // Take time step
      RK4Adaptive (t, y, h, t_err, acc, 2., rept, maxrept, hmin, hmax, 2, 0, NULL);
      Unpack (y);

      dt += h;
      
      // Update time in ms
      TIME = t * tau_A * 1.e3;

      // Output Stage6 mode locking data for IslandDynamics
      for (int j = 0; j < nres; j++)
	ww (j) = GetActualFrequency (j);
      
      for (int j = 0; j < nres; j++)
	if (ww (j) * wwo (j) < 0. && lock (j) == 0)
	  {
	    CalcRMP (t); 

	    printf ("m = %3d locks at t = %11.4e s  irmp = %11.4e kA  prmp/pi = %11.4e\n",
		    mk (j), t*tau_A, irmp, prmp /M_PI);
	 
	    fprintf (file8, "%16.9e %3d %16.9e %16.9e %16.9e %16.9e %3d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
		     q95, mk (j), rk (j), t*tau_A, irmp, prmp /M_PI, nres, q0, qa,
		     GetNaturalFrequency (j)/tau_A/1.e3, wkl (j)/tau_A/1.e3, wke (j)/tau_A/1.e3, TIME, wkn (j)/tau_A/1.e3);

	    lock (j) = 1;
	  }

      // Output Stage5 data
      if (dt > dTT)
	{
	  dt = 0.;

	  // Output reconnected fluxes
	  fprintf (file1, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file1, "%16.9e ", Psik (j));
	  fprintf (file1, "\n");

	  // Output phases of reconnected fluxes
	  fprintf (file2, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file2, "%16.9e ", atan2 (sin (phik (j) - zeta (j)), cos (phik (j) - zeta (j))) /M_PI);
	  fprintf (file2, "\n");

	  // Output island widths in r
	  fprintf (file3, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file3, "%16.9e ", GetIslandWidth (j) * R_0 /dPsiNdr (j));
	  fprintf (file3, "\n");

	  // Output ratios of island widths to linear layer widths
	  fprintf (file3a, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file3a, "%16.9e ", GetIslandWidth (j) * R_0 /dPsiNdr (j) /delk (j));
	  fprintf (file3a, "\n");

	  // Output modified natural frequencies
	  fprintf (file4, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file4, "%16.9e ", ww (j) /tau_A/1.e3);
	  fprintf (file4, "\n");

	  // Output unmodified natural frequencies
	  fprintf (file4a, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file4a, "%16.9e ", GetNaturalFrequency (j) /tau_A/1.e3);
	  fprintf (file4a, "\n");

	  // Output changes in poloidal angular velocities
	  fprintf (file4b, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file4b, "%16.9e ", GetDeltaOmegaTheta (j));
	  fprintf (file4b, "\n");

	  // Output changes in toroidal angular velocities
	  fprintf (file4c, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file4c, "%16.9e ", GetDeltaOmegaPhi (j));
	  fprintf (file4c, "\n");

	  // Output changes in Er
	  fprintf (file4d, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file4d, "%16.9e ", GetDeltaEr (j));
	  fprintf (file4d, "\n");
	  
	  // Output RMP data
	  CalcRMP (t);
	  fprintf (file5, "%16.9e %16.9e %16.9e\n", t*tau_A, irmp, prmp /M_PI);

	  // Output coil current data
	  double IU, IM, IL, PU, PM, PL;
	  CalcCoil (t, IU, IM, IL, PU, PM, PL);
	  
	  int    k       = Findk ();
	  double chikm   = chi (k);
	  double chikp   = chi (k+1);
	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));
	  double chik    = chikm * Weightm + chikp * Weightp;
	  double chik1   = chi (0);

	  fprintf (file20, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e \n", t*tau_A, IU, IM, IL, PU/M_PI, PM/M_PI, PL/M_PI, chik, chik1); 

	  // Calculate and output simulated Mirnov data
	  double sumci = 0., sumsi = 0., sumco = 0., sumso = 0.;

	  for (int j = 0; j < nres; j++)
	    {
	      sumci += epsi (j) * EEh (j, j) * chi (j) * cos (sigi (j) + zeta (j));
	      sumsi += epsi (j) * EEh (j, j) * chi (j) * sin (sigi (j) + zeta (j));
	      sumco += epso (j) * EEh (j, j) * chi (j) * cos (sigo (j) + zeta (j));
	      sumso += epso (j) * EEh (j, j) * chi (j) * sin (sigo (j) + zeta (j));
	      
	      for (int k = 0; k < nres; k++)
		{
		  sumci += epsi (j) * EEh (j, k) * Psik (k) * cos (sigi (j) + phik (k) + xih (j, k));
		  sumsi += epsi (j) * EEh (j, k) * Psik (k) * sin (sigi (j) + phik (k) + xih (j, k));
		  sumco += epso (j) * EEh (j, k) * Psik (k) * cos (sigo (j) + phik (k) + xih (j, k));
		  sumso += epso (j) * EEh (j, k) * Psik (k) * sin (sigo (j) + phik (k) + xih (j, k));
		}
	    }
	  
	  fprintf (file6, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", t*tau_A, sumci*B_0*1.e4, sumsi*B_0*1.e4, sumco*B_0*1.e4, sumso*B_0*1.e4);

	  fprintf (file9,  "%16.9e ", t*tau_A); fprintf (file10, "%16.9e ", t*tau_A); fprintf (file11, "%16.9e ", t*tau_A);
	  fprintf (file12, "%16.9e ", t*tau_A); fprintf (file13, "%16.9e ", t*tau_A); fprintf (file14, "%16.9e ", t*tau_A);
	  fprintf (file15, "%16.9e ", t*tau_A); fprintf (file16, "%16.9e ", t*tau_A); fprintf (file17, "%16.9e ", t*tau_A);
	  fprintf (file18, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    {
	      // Calculate island width in PsiN
	      double Wpk = GetIslandWidth (j);

	      // Calculate island width in r
	      double Wrk = R_0 * Wpk /dPsiNdr (j);

	      // Calculate vacuum island width in PsiN
	      double Wvk = GetVacuumIslandWidth (j);

	      // Calculate density and temperature flattening widths in r
	      double deltanek = (2./M_PI) * Wrk *Wrk*Wrk /(Wrk*Wrk + Wcrnek (j) * Wcrnek (j));
	      double deltaTek = (2./M_PI) * Wrk *Wrk*Wrk /(Wrk*Wrk + WcrTek (j) * WcrTek (j));

	      // Calculate density and temperature reductions
	      double dnek = dnedrk (j) * deltanek;
	      double dTek = dTedrk (j) * deltaTek;

	      // Output data
	      fprintf (file9,  "%16.9e ", deltanek);
	      fprintf (file10, "%16.9e ", deltaTek);
	      fprintf (file11, "%16.9e ", dnek);
	      fprintf (file12, "%16.9e ", dTek);
	      fprintf (file13, "%16.9e ", rk (j) - Wrk /2./a (j)/R_0);
	      fprintf (file14, "%16.9e ", rk (j) + Wrk /2./a (j)/R_0);
	      fprintf (file15, "%16.9e ", rk (j) - deltanek /2./a (j)/R_0);
	      fprintf (file16, "%16.9e ", rk (j) + deltanek /2./a (j)/R_0);
	      fprintf (file17, "%16.9e ", rk (j) - deltaTek /2./a (j)/R_0);
	      fprintf (file18, "%16.9e ", rk (j) + deltaTek /2./a (j)/R_0);

	      fprintf (file19, "%3d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
		       mk (j),
		       rk (j),
		       GetNaturalFrequency (j) /tau_A/1.e3,
		       GetActualFrequency (j)  /tau_A/1.e3,
		       t*tau_A,
		       (Wrk /R_0) /a (j),
		       PsiN (j),
		       Wpk,
		       Wvk,
		       deltanek * dPsiNdr (j) /R_0, deltaTek * dPsiNdr (j) /R_0);
	    }
	  fprintf (file9,  "\n"); fprintf (file10, "\n"); fprintf (file11, "\n"); fprintf (file12, "\n");
	  fprintf (file13, "\n"); fprintf (file14, "\n"); fprintf (file15, "\n"); fprintf (file16, "\n");
	  fprintf (file17, "\n"); fprintf (file18, "\n");

	  if (cnt%10 == 0)
	    printf ("t(ms) = %11.4e  h(ms) = %11.4e  h/tau_A = %11.4e  irmp(kA) = %11.4e  prmp/pi = %11.4e\n", t*tau_A*1.e3, h*tau_A*1.e3, h, irmp, prmp /M_PI);
	  cnt++;
	  
	  fflush (file1);  fflush (file2);  fflush (file3);  fflush (file4);  fflush (file5);
	  fflush (file6);  fflush (file8);  fflush (file9);  fflush (file10); fflush (file11);
	  fflush (file12); fflush (file13); fflush (file14); fflush (file15); fflush (file16);
	  fflush (file17); fflush (file18); fflush (file3a); fflush (file4a); fflush (file19);
	  fflush (file20); fflush (file4b); fflush (file4c); fflush (file4d); 
	}
    }
  while (t < Tend);

  printf ("t(ms) = %11.4e  h(ms) = %11.4e  h/tau_A = %11.4e  irmp(kA) = %11.4e  prmp/pi = %11.4e\n", t*tau_A*1.e3, h*tau_A*1.e3, h, irmp, prmp /M_PI);

  fclose (file1);  fclose (file2);  fclose (file3);  fclose (file4);  fclose (file5);
  fclose (file6);  fclose (file8);  fclose (file9);  fclose (file10); fclose (file11);
  fclose (file12); fclose (file13); fclose (file14); fclose (file15); fclose (file16);
  fclose (file17); fclose (file18); fclose (file3a); fclose (file4a); fclose (file19);
  fclose (file20); fclose (file4b); fclose (file4c); fclose (file4d); 

  // Save calculation for restart
  Save ();

  // Output Stage6 data
  FILE* filex = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/omega0.txt");
  for (int j = 0; j < nres; j++)
    fprintf (filex, "%3d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mk (j), rk (j), wkl (j) /tau_A/1.e3, wke (j) /tau_A/1.e3, wkn (j) /tau_A/1.e3, GetNaturalFrequency (j) /tau_A/1.e3, TIME, q95);
  fclose (filex);

  FILE* filew = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/omega.txt");
  for (int j = 0; j < nres; j++)
    {
      // Calculate island width in PsiN
      double Wpk = GetIslandWidth (j);
       
      // Calculate island width in r
      double Wrk = R_0 * Wpk /dPsiNdr (j);

      // Calculate vacuum island width in PsiN
      double Wvk = GetVacuumIslandWidth (j);

      // Calculate density and temperature flattening widths in PsiN
      double deltanek = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + Wcrnek (j) * Wcrnek (j));
      double deltaTek = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + WcrTek (j) * WcrTek (j));

      // Calculate radial magnetic field at edge
      double br_unrc = 1.e4 * EEh (j, j) * chi (j) * (pow (rk (j), - double (mk (j))) - pow (rk (j), + double (mk (j)))) * fabs (B_0) /2./a (j);
      double br_full = 1.e4 * EEh (j, j) * chi (j) *  pow (rk (j), - double (mk (j)))                                    * fabs (B_0) /2./a (j);

      // Calculate changes in poloidal angular velocity, toroidal angular velocity, toroidal velocity, ExB velocity, and radial electric field
      double Omegat = GetDeltaOmegaTheta (j);
      double Omegap = GetDeltaOmegaPhi   (j);
      double Vp     = GetDeltaVPhi       (j);
      double VEB    = GetDeltaVEB        (j);
      double Er     = GetDeltaEr         (j);
      double omega  = GetDeltaOmega      (j);

      fprintf (filew, "%3d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	       mk (j),
	       rk (j),
	       GetNaturalFrequency (j) /tau_A/1.e3,
	       GetActualFrequency  (j) /tau_A/1.e3,
	       TIME,
	       (Wrk /R_0) /a (j),
	       PsiN (j),
	       Wpk,
	       Wvk,
	       deltanek, deltaTek, q95, br_unrc, br_full,
	       Omegat, Omegap, Vp, VEB, Er, omega);
    }
  fclose (filew);

  FILE* filewx = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/deltap.txt");
  double deltap = 0.;
  for (int j = 0; j < nres; j++)
    {
      // Calculate island width in PsiN
      double Wpk = GetIslandWidth (j);
       
      // Calculate island width in r
      double Wrk = R_0 * Wpk /dPsiNdr (j);

      // Calculate density and temperature flattening widths in PsiN
      double deltanek = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + Wcrnek (j) * Wcrnek (j));
      double deltaTek = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + WcrTek (j) * WcrTek (j));
      double deltaTik = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + WcrTik (j) * WcrTik (j));

      // Calculate pressure decrement
      double deltapk;
      if (HIGH)
	{
	  deltapk =
	      deltanek * Factor1 (j) + deltanek * Wpk*Wpk * Factor5 (j) + deltanek*deltanek*deltanek * Factor9  (j)
	    + deltaTek * Factor2 (j) + deltaTek * Wpk*Wpk * Factor6 (j) + deltaTek*deltaTek*deltaTek * Factor10 (j)
	    + deltanek * Factor3 (j) + deltanek * Wpk*Wpk * Factor7 (j) + deltanek*deltanek*deltanek * Factor11 (j)
	    + deltaTik * Factor4 (j) + deltaTik * Wpk*Wpk * Factor8 (j) + deltaTik*deltaTik*deltaTik * Factor12 (j);
	}
      else
	{
	  deltapk =
	      deltanek * Factor1 (j)
	    + deltaTek * Factor2 (j)
	    + deltanek * Factor3 (j)
	    + deltaTik * Factor4 (j);
	}

      // Calculate cumulative pressure decrement
      deltap += deltapk;
    }

  fprintf (filewx, "%16.9e %16.9e %16.9e %16.9e\n", TIME, q95, deltap /P0/Pped, deltap /P0);

  fclose (filewx);

  // Output coil optimization data
  double IU, IM, IL, PU, PM, PL;
  CalcCoil (t, IU, IM, IL, PU, PM, PL);

  int    k       = Findk ();
  double chikm   = chi (k);
  double chikp   = chi (k+1);
  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));
  double chik    = chikm * Weightm + chikp * Weightp;
  double chik1   = chi (0);

  FILE* fileco = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/opt.txt");
  fprintf (fileco, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", TIME, q95, IU, IM, IL, PU/M_PI, PM/M_PI, PL/M_PI, chik, chik1); 
  fclose (fileco);

  // Output Mirnov data
  FILE* fileww = OpenFilea ((char*) "../IslandDynamics/Outputs/Stage6/mirnov.txt");
  double sumci = 0., sumsi = 0., sumco = 0., sumso = 0.;

  for (int j = 0; j < nres; j++)
    {
      sumci += epsi (j) * EEh (j, j) * chi (j) * cos (sigi (j) + zeta (j));
      sumsi += epsi (j) * EEh (j, j) * chi (j) * sin (sigi (j) + zeta (j));
      sumco += epso (j) * EEh (j, j) * chi (j) * cos (sigo (j) + zeta (j));
      sumso += epso (j) * EEh (j, j) * chi (j) * sin (sigo (j) + zeta (j));
      
      for (int k = 0; k < nres; k++)
	{
	  sumci += epsi (j) * EEh (j, k) * Psik (k) * cos (sigi (j) + phik (k) + xih (j, k));
	  sumsi += epsi (j) * EEh (j, k) * Psik (k) * sin (sigi (j) + phik (k) + xih (j, k));
	  sumco += epso (j) * EEh (j, k) * Psik (k) * cos (sigo (j) + phik (k) + xih (j, k));
	  sumso += epso (j) * EEh (j, k) * Psik (k) * sin (sigo (j) + phik (k) + xih (j, k));
	}
    }
  
  fprintf (fileww, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", t*tau_A*1.e3, sumci*B_0*1.e4, sumsi*B_0*1.e4, sumco*B_0*1.e4, sumso*B_0*1.e4);
  fclose (fileww);
  
  // Output magnetic island chains versus theta
  FILE* filep = OpenFilew ((char*) "Outputs/Stage5/islandt.txt");
  int IMAX = 2000;
  for (int i = 0; i < IMAX; i++)
    {
      double theta = 2.*M_PI * double (i) /double (IMAX - 1);

      double Xminus, Xplus;
      for (int j = 0; j < nres; j++)
	{
	  GetIslandLimits (j, Psik (j) * cos (double (mk (j)) * theta - phik (j)), Xminus, Xplus);
	  fprintf (filep, "%d %e %e %e\n", mk (j), theta/M_PI, PsiN (j) + Xminus, PsiN (j) + Xplus);
	}
    }
  fclose (filep);

  // Output magnetic island chains versus phi
  FILE* fileq = OpenFilew ((char*) "Outputs/Stage5/islandp.txt");
  for (int i = 0; i < IMAX; i++)
    {
      double phi = 2.*M_PI * double (i) /double (IMAX - 1);

      double Xminus, Xplus;
      for (int j = 0; j < nres; j++)
	{
	  GetIslandLimits (j, Psik (j) * cos (double (mk (j)) * M_PI - double (ntor (j)) * phi - phik (j)), Xminus, Xplus);
	  fprintf (fileq, "%d %e %e %e\n", mk (j), phi/M_PI, PsiN (j) + Xminus, PsiN (j) + Xplus);
	}
    }
  fclose (fileq);

  // Output vacuum magnetic island chains versus theta
  FILE* filepv = OpenFilew ((char*) "Outputs/Stage5/islandtv.txt");
  for (int i = 0; i < IMAX; i++)
    {
      double theta = 2.*M_PI * double (i) /double (IMAX - 1);

      double Xminus, Xplus;
      for (int j = 0; j < nres; j++)
	{
	  GetIslandLimits (j, chi (j) * cos (double (mk (j)) * theta - zeta (j)), Xminus, Xplus);
	  fprintf (filepv, "%d %e %e %e\n", mk (j), theta/M_PI, PsiN (j) + Xminus, PsiN (j) + Xplus);
	}
    }
  fclose (filepv);

  // Output vacuum magnetic island chains versus phi
  FILE* fileqv = OpenFilew ((char*) "Outputs/Stage5/islandpv.txt");
  for (int i = 0; i < IMAX; i++)
    {
      double phi = 2.*M_PI * double (i) /double (IMAX - 1);

      double Xminus, Xplus;
      for (int j = 0; j < nres; j++)
	{
	  GetIslandLimits (j, chi (j) * cos (double (mk (j)) * M_PI - double (ntor (j)) * phi - zeta (j)), Xminus, Xplus);
	  fprintf (fileqv, "%d %e %e %e\n", mk (j), phi/M_PI, PsiN (j) + Xminus, PsiN (j) + Xplus);
	}
    }
  fclose (fileqv);
}

// ###################################
// Function to calculate irmp and prmp
// ###################################
void Phase::CalcRMP (double t)
{
  if (t <= TT (0))
    {
      irmp = ICTRL[0];
      prmp = PCTRL[0];
    }
  else if (t >= TT (NCTRL-1))
    {
      irmp = ICTRL[NCTRL-1];
      prmp = PCTRL[NCTRL-1];
    }
  else
    {
      for (int i = 0; i < NCTRL-1; i++)
	{
	  if ((t > TT (i)) && (t < TT (i+1)))
	    {
	      if (IFLA)
		irmp = IRMP;
	      else
		irmp = ((t - TT (i)) * ICTRL[i+1] + (TT (i+1) - t) * ICTRL[i]) /(TT (i+1) - TT (i));
	      prmp = ((t - TT (i)) * PCTRL[i+1] + (TT (i+1) - t) * PCTRL[i]) /(TT (i+1) - TT (i));
	    }
	}
    }
}

// ##############################################
// Function to calculate coil currents and phases
// ##############################################
void Phase::CalcCoil (double t, double& IU, double& IM, double& IL, double& PU, double& PM, double& PL)
{
  CalcRMP (t);

  if (COPT == 0)
    {
      if (MID == 3)
	{
	  IU = irmp;
	  IM = irmp;
	  IL = irmp;
	  PU = - prmp;
	  PM = 0.;
	  PL = + prmp;
	}
      else if (MID == 2)
	{
	  IL = irmp;
	  IM = 0.;
	  IU = irmp;
	  PU = - prmp /2.;
	  PM = 0.;
	  PL = + prmp /2.;
	}
      else
	{
	  IU = 0.;
	  IM = 0.;
	  IL = irmp;
	  PU = 0.;
	  PM = 0.;
	  PL = 0.;
	}
    }
  else if (COPT == 1)
    {
      if (MID == 3)
	{
	  int k = Findk ();
	
	  gsl_complex chikmL = gsl_vector_complex_get (ChiL, k);
	  gsl_complex chikmU = gsl_vector_complex_get (ChiU, k);
	  gsl_complex chikmM = gsl_vector_complex_get (ChiM, k);

	  gsl_complex chikpL = gsl_vector_complex_get (ChiL, k+1);
	  gsl_complex chikpU = gsl_vector_complex_get (ChiU, k+1);
	  gsl_complex chikpM = gsl_vector_complex_get (ChiM, k+1);

	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));

	  chikmL = gsl_complex_mul_real (chikmL, Weightm);
	  chikmU = gsl_complex_mul_real (chikmU, Weightm);
	  chikmM = gsl_complex_mul_real (chikmM, Weightm);
					 
	  chikpL = gsl_complex_mul_real (chikpL, Weightp);
	  chikpU = gsl_complex_mul_real (chikpU, Weightp);
	  chikpM = gsl_complex_mul_real (chikpM, Weightp);

	  gsl_complex chikL = gsl_complex_add (chikmL, chikpL);
	  gsl_complex chikU = gsl_complex_add (chikmU, chikpU);
	  gsl_complex chikM = gsl_complex_add (chikmM, chikpM);

	  gsl_complex chikLa = gsl_complex_conjugate (chikL);
	  gsl_complex chikMa = gsl_complex_conjugate (chikM);

	  gsl_complex xkUM = gsl_complex_mul (chikU, chikMa);
	  gsl_complex xkUL = gsl_complex_mul (chikU, chikLa);
	  gsl_complex xkML = gsl_complex_mul (chikM, chikLa);

	  double XkUM = gsl_complex_abs (xkUM);
	  double XkUL = gsl_complex_abs (xkUL);
	  double XkML = gsl_complex_abs (xkML);
	  
	  double gammakUM = gsl_complex_arg (xkUM);
	  double gammakUL = gsl_complex_arg (xkUL);
	  double gammakML = gsl_complex_arg (xkML);
	  
	  double Delta = FindMax (XkUM, XkUL, XkML, gammakUM, gammakUL, gammakML);
	  
	  IU = irmp;
	  IM = irmp;
	  IL = irmp;
  	  PU =  - Delta;
	  PM =  0.;
	  PL =  + Delta;
	}
      else if (MID == 2)
	{
	  int k = Findk ();
	
	  gsl_complex chikmL = gsl_vector_complex_get (ChiL, k);
	  gsl_complex chikmU = gsl_vector_complex_get (ChiU, k);

	  gsl_complex chikpL = gsl_vector_complex_get (ChiL, k+1);
	  gsl_complex chikpU = gsl_vector_complex_get (ChiU, k+1);

	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));

	  chikmL = gsl_complex_mul_real (chikmL, Weightm);
	  chikmU = gsl_complex_mul_real (chikmU, Weightm);

	  chikpL = gsl_complex_mul_real (chikpL, Weightp);
	  chikpU = gsl_complex_mul_real (chikpU, Weightp);

	  gsl_complex chikL = gsl_complex_add (chikmL, chikpL);
	  gsl_complex chikU = gsl_complex_add (chikmU, chikpU);

	  gsl_complex chikLa = gsl_complex_conjugate (chikL);

	  gsl_complex xkUL = gsl_complex_mul (chikU, chikLa);

	  double gammakUL = gsl_complex_arg (xkUL);

	  double Delta = gammakUL;

  	  IU = irmp;
	  IM = 0.;
	  IL = irmp;
	  PU =  - Delta /2.;
	  PM =  0.;
	  PL =  + Delta /2.;
	}
      else 
	{
	  IU = 0.;
	  IM = 0.;
	  IL = irmp;
  	  PU = 0.;
	  PM = 0.;
	  PL = 0.;
	}
    }
  else if (COPT == 2)
    {
      if (MID == 3)
	{
	  int k = Findk ();
	
	  gsl_complex chikmL = gsl_vector_complex_get (ChiL, k);
	  gsl_complex chikmU = gsl_vector_complex_get (ChiU, k);
	  gsl_complex chikmM = gsl_vector_complex_get (ChiM, k);

	  gsl_complex chikpL = gsl_vector_complex_get (ChiL, k+1);
	  gsl_complex chikpU = gsl_vector_complex_get (ChiU, k+1);
	  gsl_complex chikpM = gsl_vector_complex_get (ChiM, k+1);

	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));

	  chikmL = gsl_complex_mul_real (chikmL, Weightm);
	  chikmU = gsl_complex_mul_real (chikmU, Weightm);
	  chikmM = gsl_complex_mul_real (chikmM, Weightm);
					 
	  chikpL = gsl_complex_mul_real (chikpL, Weightp);
	  chikpU = gsl_complex_mul_real (chikpU, Weightp);
	  chikpM = gsl_complex_mul_real (chikpM, Weightp);

	  gsl_complex chikL = gsl_complex_add (chikmL, chikpL);
	  gsl_complex chikU = gsl_complex_add (chikmU, chikpU);
	  gsl_complex chikM = gsl_complex_add (chikmM, chikpM);

	  gsl_complex iL = gsl_complex_conjugate (chikL);
	  gsl_complex iU = gsl_complex_conjugate (chikU);
	  gsl_complex iM = gsl_complex_conjugate (chikM);

	  double      argM = gsl_complex_arg (iM);
 	  gsl_complex eiM  = gsl_complex_polar (1., argM);

	  iL = gsl_complex_div (iL, eiM);
	  iU = gsl_complex_div (iU, eiM);
	  iM = gsl_complex_div (iM, eiM);

	  double I = (gsl_complex_abs (iL) + gsl_complex_abs (iU) + gsl_complex_abs (iM)) /3.;

	  IU = irmp * gsl_complex_abs (iU) /I;
	  IM = irmp * gsl_complex_abs (iM) /I;
	  IL = irmp * gsl_complex_abs (iL) /I;
	  PU = gsl_complex_arg (iU);
	  PM = gsl_complex_arg (iM);
	  PL = gsl_complex_arg (iL);
	}
      else if (MID == 2)
	{
	  int k = Findk ();
	
	  gsl_complex chikmL = gsl_vector_complex_get (ChiL, k);
	  gsl_complex chikmU = gsl_vector_complex_get (ChiU, k);

	  gsl_complex chikpL = gsl_vector_complex_get (ChiL, k+1);
	  gsl_complex chikpU = gsl_vector_complex_get (ChiU, k+1);

	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));

	  chikmL = gsl_complex_mul_real (chikmL, Weightm);
	  chikmU = gsl_complex_mul_real (chikmU, Weightm);
					 
	  chikpL = gsl_complex_mul_real (chikpL, Weightp);
	  chikpU = gsl_complex_mul_real (chikpU, Weightp);

	  gsl_complex chikL = gsl_complex_add (chikmL, chikpL);
	  gsl_complex chikU = gsl_complex_add (chikmU, chikpU);

	  gsl_complex iL = gsl_complex_conjugate (chikL);
	  gsl_complex iU = gsl_complex_conjugate (chikU);

	  double      argM = (gsl_complex_arg (iL) + gsl_complex_arg (iU)) /2.;
	  gsl_complex eiM  = gsl_complex_polar (1., argM);

	  iL = gsl_complex_div (iL, eiM);
	  iU = gsl_complex_div (iU, eiM);

	  double I = (gsl_complex_abs (iL) + gsl_complex_abs (iU)) /2.;

	  IU = irmp * gsl_complex_abs (iU) /I;
	  IM = 0.;
	  IL = irmp * gsl_complex_abs (iL) /I;
	  PU = gsl_complex_arg (iU);
	  PM = 0.;
	  PL = gsl_complex_arg (iL);
	}
      else 
	{
	  IU = 0.;
	  IM = 0.;
	  IL = irmp;
  	  PU = 0.;
	  PM = 0.;
	  PL = 0.;
	}
    }
  else if (COPT == 3)
    {
      if (MID == 3)
	{
	  int k = Findk ();
	
	  gsl_complex chikmL = gsl_vector_complex_get (ChiL, k);
	  gsl_complex chikmU = gsl_vector_complex_get (ChiU, k);
	  gsl_complex chikmM = gsl_vector_complex_get (ChiM, k);

	  gsl_complex chikpL = gsl_vector_complex_get (ChiL, k+1);
	  gsl_complex chikpU = gsl_vector_complex_get (ChiU, k+1);
	  gsl_complex chikpM = gsl_vector_complex_get (ChiM, k+1);

	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));

	  chikmL = gsl_complex_mul_real (chikmL, Weightm);
	  chikmU = gsl_complex_mul_real (chikmU, Weightm);
	  chikmM = gsl_complex_mul_real (chikmM, Weightm);
					 
	  chikpL = gsl_complex_mul_real (chikpL, Weightp);
	  chikpU = gsl_complex_mul_real (chikpU, Weightp);
	  chikpM = gsl_complex_mul_real (chikpM, Weightp);

	  gsl_complex chikL = gsl_complex_add (chikmL, chikpL);
	  gsl_complex chikU = gsl_complex_add (chikmU, chikpU);
	  gsl_complex chikM = gsl_complex_add (chikmM, chikpM);
	 
	  gsl_complex chi1L = gsl_vector_complex_get (ChiL, 0);
	  gsl_complex chi1U = gsl_vector_complex_get (ChiU, 0);
	  gsl_complex chi1M = gsl_vector_complex_get (ChiM, 0);

	  gsl_complex lambda = gsl_complex_mul (gsl_complex_conjugate (chikL), chi1L);
	  lambda             = gsl_complex_add (lambda, gsl_complex_mul (gsl_complex_conjugate (chikU), chi1U));
	  lambda             = gsl_complex_add (lambda, gsl_complex_mul (gsl_complex_conjugate (chikM), chi1M));
	  lambda             = gsl_complex_div_real (lambda, gsl_complex_abs2 (chi1L) + gsl_complex_abs2 (chi1U) + gsl_complex_abs2 (chi1M));
	  lambda             = gsl_complex_mul_real (lambda, CORE);
	  
	  gsl_complex iL = gsl_complex_sub (gsl_complex_conjugate (chikL), gsl_complex_mul (lambda, gsl_complex_conjugate (chi1L)));
	  gsl_complex iU = gsl_complex_sub (gsl_complex_conjugate (chikU), gsl_complex_mul (lambda, gsl_complex_conjugate (chi1U)));
	  gsl_complex iM = gsl_complex_sub (gsl_complex_conjugate (chikM), gsl_complex_mul (lambda, gsl_complex_conjugate (chi1M)));

	  double      argM = gsl_complex_arg (iM);
 	  gsl_complex eiM  = gsl_complex_polar (1., argM);

	  iL = gsl_complex_div (iL, eiM);
	  iU = gsl_complex_div (iU, eiM);
	  iM = gsl_complex_div (iM, eiM);

	  double I = (gsl_complex_abs (iL) + gsl_complex_abs (iU) + gsl_complex_abs (iM)) /3.;

	  IU = irmp * gsl_complex_abs (iU) /I;
	  IM = irmp * gsl_complex_abs (iM) /I;
	  IL = irmp * gsl_complex_abs (iL) /I;
	  PU = gsl_complex_arg (iU);
	  PM = gsl_complex_arg (iM);
	  PL = gsl_complex_arg (iL);
	}
      else if (MID == 2)
	{
	  int k = Findk ();
	
	  gsl_complex chikmL = gsl_vector_complex_get (ChiL, k);
	  gsl_complex chikmU = gsl_vector_complex_get (ChiU, k);

	  gsl_complex chikpL = gsl_vector_complex_get (ChiL, k+1);
	  gsl_complex chikpU = gsl_vector_complex_get (ChiU, k+1);

	  double Weightm = (PsiN (k+1) - PSIPED)   /(PsiN (k+1) - PsiN (k));
	  double Weightp = (PSIPED     - PsiN (k)) /(PsiN (k+1) - PsiN (k));

	  chikmL = gsl_complex_mul_real (chikmL, Weightm);
	  chikmU = gsl_complex_mul_real (chikmU, Weightm);
					 
	  chikpL = gsl_complex_mul_real (chikpL, Weightp);
	  chikpU = gsl_complex_mul_real (chikpU, Weightp);

	  gsl_complex chikL = gsl_complex_add (chikmL, chikpL);
	  gsl_complex chikU = gsl_complex_add (chikmU, chikpU);

	  gsl_complex chi1L = gsl_vector_complex_get (ChiL, 0);
	  gsl_complex chi1U = gsl_vector_complex_get (ChiU, 0);

	  gsl_complex lambda = gsl_complex_mul (gsl_complex_conjugate (chikL), chi1L);
	  lambda             = gsl_complex_add (lambda, gsl_complex_mul (gsl_complex_conjugate (chikU), chi1U));
	  lambda             = gsl_complex_div_real (lambda, gsl_complex_abs2 (chi1L) + gsl_complex_abs2 (chi1U));
	  lambda             = gsl_complex_mul_real (lambda, CORE);
	  
	  gsl_complex iL = gsl_complex_sub (gsl_complex_conjugate (chikL), gsl_complex_mul (lambda, gsl_complex_conjugate (chi1L)));
	  gsl_complex iU = gsl_complex_sub (gsl_complex_conjugate (chikU), gsl_complex_mul (lambda, gsl_complex_conjugate (chi1U)));

	  double      argM = (gsl_complex_arg (iL) + gsl_complex_arg (iU)) /2.;
	  gsl_complex eiM  = gsl_complex_polar (1., argM);

	  iL = gsl_complex_div (iL, eiM);
	  iU = gsl_complex_div (iU, eiM);

	  double I = (gsl_complex_abs (iL) + gsl_complex_abs (iU)) /2.;

	  IU = irmp * gsl_complex_abs (iU) /I;
	  IM = 0.;
	  IL = irmp * gsl_complex_abs (iL) /I;
	  PU = gsl_complex_arg (iU);
	  PM = 0.;
	  PL = gsl_complex_arg (iL);
	}
      else 
	{
	  IU = 0.;
	  IM = 0.;
	  IL = irmp;
  	  PU = 0.;
	  PM = 0.;
	  PL = 0.;
	}
    }
  else
    {
      printf ("PHASE::Error unknown option COPT = %2d\n");
      exit (1);
    }
}

// ################################################################
// Function to find resonant surfaces that straddle top of pedestal
// ################################################################
int Phase::Findk ()
{
  int k = -1;

  for (int j = 0; j < nres-1; j++)
    if ((PsiN (j) - PSIPED) * (PsiN (j+1) - PSIPED) <= 0.)
      k = j;

  if (k == -1)
    k = nres - 2;

  return k;
}

// ##########################################################
// Function to find maximum of restricted three-coil function
// ##########################################################
double Phase::FindMax (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML)
{
  double Delta = 0., fun, deriv, dderiv, max = -1.e6;
  
  for (int i = 0; i < 360; i++)
    {
      double delta = (double (i) /360.) * 2.*M_PI;

      ThreeCoil (XUM, XUL, XML, gammaUM, gammaUL, gammaML, delta, fun, deriv, dderiv);

      if (fun > max)
	{
	  Delta = delta;
	  max   = fun;
	}
    }

  for (int i = 0; i < 5; i++)
    {
      ThreeCoil (XUM, XUL, XML, gammaUM, gammaUL, gammaML, Delta, fun, deriv, dderiv);

      Delta -= deriv /dderiv; 
    }

  return Delta;
}

// ##########################################################
// Function to find minimum of restricted three-coil function
// ##########################################################
double Phase::FindMin (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML)
{
  double Delta = 0., fun, deriv, dderiv, min = 1.e6;
  
  for (int i = 0; i < 360; i++)
    {
      double delta = (double (i) /360.) * 2.*M_PI;

      ThreeCoil (XUM, XUL, XML, gammaUM, gammaUL, gammaML, delta, fun, deriv, dderiv);

      if (fun < min)
	{
	  Delta = delta;
	  min   = fun;
	}
    }

  for (int i = 0; i < 5; i++)
    {
      ThreeCoil (XUM, XUL, XML, gammaUM, gammaUL, gammaML, Delta, fun, deriv, dderiv);

      Delta -= deriv /dderiv; 
    }

  return Delta;
}

// #######################################################################
// Function to evaluate restricted three-coil function and its derivatives
// #######################################################################
void Phase::ThreeCoil (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML, double Delta, double& fun, double& deriv, double& dderiv)
{
  fun    =   XUM * cos (Delta - gammaUM) +      XUL * cos (2.*Delta - gammaUL) + XML * cos (Delta - gammaML);
  deriv  = - XUM * sin (Delta - gammaUM) - 2. * XUL * sin (2.*Delta - gammaUL) - XML * sin (Delta - gammaML);
  dderiv = - XUM * cos (Delta - gammaUM) - 4. * XUL * cos (2.*Delta - gammaUL) - XML * cos (Delta - gammaML);
}

// ##########################################
// Function to calculate chi and zeta vectors
// ##########################################
void Phase::CalcChiZeta (double t)
{
  // Get coil currents and phases
  double IU, IM, IL, PU, PM, PL;
  CalcCoil (t, IU, IM, IL, PU, PM, PL);

  // Calculate complex coil currents
  double      one = 1.;
  gsl_complex eiL = gsl_complex_polar    (one, PL);
  gsl_complex eiU = gsl_complex_polar    (one, PU);
  gsl_complex eiM = gsl_complex_polar    (one, PM);
  gsl_complex IIL = gsl_complex_mul_real (eiL, IL);
  gsl_complex IIU = gsl_complex_mul_real (eiU, IU);
  gsl_complex IIM = gsl_complex_mul_real (eiM, IM);
  gsl_complex h;
 
  // Calculate chi and zeta values
  for (int j = 0; j < nres; j++)
    {
      if (MID == 3)
	{
	  gsl_complex hl = gsl_vector_complex_get (ChiL, j);
	  gsl_complex hu = gsl_vector_complex_get (ChiU, j);
	  gsl_complex hm = gsl_vector_complex_get (ChiM, j);

	  hl = gsl_complex_mul (IIL, hl);
	  hu = gsl_complex_mul (IIU, hu);
	  hm = gsl_complex_mul (IIM, hm);

	  h = gsl_complex_add (hl, hu);
	  h = gsl_complex_add (h,  hm);
	}
      else if (MID == 2)
	{
	  gsl_complex hl = gsl_vector_complex_get (ChiL, j);
	  gsl_complex hu = gsl_vector_complex_get (ChiU, j);

	  hl = gsl_complex_mul (IIL, hl);
	  hu = gsl_complex_mul (IIU, hu);

	  h = gsl_complex_add (hl, hu);
	}
      else
	{
	  gsl_complex hl = gsl_vector_complex_get (ChiL, j);

	  hl = gsl_complex_mul (IIL, hl);

	  h = hl;
	}
    
      chi  (j) =   gsl_complex_abs (h);
      zeta (j) = - gsl_complex_arg (h);

      // Limit maximum Chirikov parameter for vacuum islands
      double Wvac = GetVacuumIslandWidth (j);
      chi (j)     = Wvac*Wvac /16. /A1 (j);
    }
}

// ########################################################
// Function to pack simulation variables into single vector
// ########################################################
void Phase::Pack (Array<double,1> y)
{
  int cnt = 0;
  for (int j = 0; j < nres; j++)
    {
      y (cnt) = Xk (j); cnt++;
      y (cnt) = Yk (j); cnt++;
      for (int i = 0; i < NFLOW; i++)
	{
	  y (cnt) = alphakp (j, i); cnt++;
	  y (cnt) = betakp  (j, i); cnt++;
	}
    }
}

// ##########################################################
// Function to unpack simulation variables from single vector
// ##########################################################
void Phase::Unpack (Array<double,1> y)
{
  int cnt = 0;
  for (int j = 0; j < nres; j++)
    {
      Xk   (j) = y (cnt); cnt++;
      Yk   (j) = y (cnt); cnt++;
      Psik (j) = sqrt (Xk (j) * Xk (j) + Yk (j) * Yk (j));
      phik (j) = atan2 (Yk (j), Xk (j));
      for (int i = 0; i < NFLOW; i++)
	{
	  alphakp (j, i) = y (cnt); cnt++;
	  betakp  (j, i) = y (cnt); cnt++;
	}
    }
}

// ####################################################
// Function to pack right-hand sides into single vector
// ####################################################
void Phase::PackRhs (Array<double,1> XkRHS, Array<double,1> YkRHS,
		     Array<double,2> alphakpRHS, Array<double,2> betakpRHS,
		     Array<double,1> dydt)
{
  int cnt = 0;
  for (int j = 0; j < nres; j++)
    {
      dydt (cnt) = XkRHS (j); cnt++;
      dydt (cnt) = YkRHS (j); cnt++;
      for (int i = 0; i < NFLOW; i++)  
	{
	  dydt (cnt) = alphakpRHS (j, i); cnt++;
	  dydt (cnt) = betakpRHS  (j, i); cnt++;
	}
    }

  for (int j = 0; j < cnt-1; j++)
    {
      if (isnan (dydt(j)))
	{
	  printf ("PHASE: Error dydt(%4d) = NaN\n", j);
	  exit (1);
	}
    }
}

// #######################################
// Function to calculate natural frequency
// #######################################
double Phase::GetNaturalFrequency (int j)
{
  if (LIN)
    {
      double om = wkl (j);

      if (om > omegamax*tau_A*1.e3)
	om =   omegamax*tau_A*1.e3;
      else if (om < - omegamax*tau_A*1.e3)
	om = - omegamax*tau_A*1.e3;

      return om;
    }
  else
    {
      double om;
      
      if (FREQ == 0)
	{
	  double w = (0.8227/2.) * 4. * R_0 * fack (j) * sqrt (fabs (Psik (j))) /delk (j);
	  
	  om = (wkl (j) + wkn (j) * w) /(1. + w);
	}
      else if (FREQ == 1)
	{
	  double w = (0.8227/2.) * 4. * R_0 * fack (j) * sqrt (fabs (Psik (j))) /delk (j);

	  om = (wkl (j) + (wke (j) - wkl (j) - wkn (j)) * w + wkn (j) * w*w) /(1. - w + w*w);
	}
      else
	om = FFAC * wkl (j) +  (1. - FFAC) * wke (j);

      if (om > omegamax*tau_A*1.e3)
	om =   omegamax*tau_A*1.e3;
      else if (om < - omegamax*tau_A*1.e3)
	om = - omegamax*tau_A*1.e3;
      
      return om;
    }
}

// ######################################
// Function to calculate actual frequency
// ######################################
double Phase::GetActualFrequency (int j)
{
  double sum = GetNaturalFrequency (j);

  for (int k = 0; k < nres; k++)
    for (int i = 0; i < NFLOW; i++)
      sum -= alphakp (k, i) * natp (j, k, i) + betakp (k, i) * natt (j, k, i);

  return sum;
}

// ##########################################################
// Function to calculate change in natural frequency (krad/s)
// ##########################################################
double Phase::GetDeltaOmega (int j)
{
  double sum = 0.;

  for (int k = 0; k < nres; k++)
    for (int i = 0; i < NFLOW; i++)
      sum -= alphakp (k, i) * natp (j, k, i) + betakp (k, i) * natt (j, k, i);

  return sum /tau_A/1.e3;
}

// ##################################################################
// Function to calculate change in poloidal angular velocity (krad/s)
// ##################################################################
double Phase::GetDeltaOmegaTheta (int j)
{
  double sum = 0.;

  for (int k = 0; k < nres; k++)
    for (int i = 0; i < NFLOW; i++)
      sum += alphakp (k, i) * natp (j, k, i);

  return - sum /double (mk (j)) /tau_A/1.e3;
}

// ##################################################################
// Function to calculate change in toroidal angular velocity (krad/s)
// ##################################################################
double Phase::GetDeltaOmegaPhi (int j)
{
  double sum = 0.;

  for (int k = 0; k < nres; k++)
    for (int i = 0; i < NFLOW; i++)
      sum += betakp (k, i) * natt (j, k, i);

  return sum /double (ntor (j)) /tau_A/1.e3;
}

// ########################################################
// Function to calculate change in toroidal velocity (km/s)
// ########################################################
double Phase::GetDeltaVPhi (int j)
{
  double value = GetDeltaOmegaPhi (j);

  return R_0 * value;
}

// ########################################################
// Function to calculate change in parallel velocity (km/s)
// ########################################################
double Phase::GetDeltaVParallel (int j)
{
  double delta_theta    = GetDeltaOmegaTheta (j);
  double delta_phi      = GetDeltaOmegaPhi   (j);
  double delta_parallel = (C2k (j) * delta_phi + (1. - C2k (j)) * qk (j) * delta_theta) /C1k (j);

  return R_0 * delta_parallel /gk (j);
}

// ###################################################
// Function to calculate change in ExB velocity (km/s)
// ###################################################
double Phase::GetDeltaVEB (int j)
{
  double delta_theta = GetDeltaOmegaTheta (j);
  double delta_phi   = GetDeltaOmegaPhi   (j);
  double delta_EB    = delta_phi - qk (j) * delta_theta;

  return - R_0 * C2k (j) * a (j) * rk (j) * delta_EB /qk (j);
}

// #########################################
// Function to calculate change in Er (kV/m)
// #########################################
double Phase::GetDeltaEr (int j)
{
  double delta_theta = GetDeltaOmegaTheta (j);
  double delta_phi   = GetDeltaOmegaPhi   (j);
  double delta_EB    = delta_phi - qk (j) * delta_theta;

  return R_0 * B_0 * dPsiNdr (j) * delta_EB;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Phase::Rhs (double t, Array<double,1>& y, Array<double,1>& dydt)
{
  CalcChiZeta (t);
  Unpack      (y);

  Array<double,1> Cosk (nres);
  Array<double,1> Sink (nres);
  Array<double,1> sink (nres);

  for (int j = 0; j < nres; j++)
    {
      Cosk (j) = EEh (j, j) * chi (j) * cos (zeta (j));
      Sink (j) = EEh (j, j) * chi (j) * sin (zeta (j));
      sink (j) = EEh (j, j) * chi (j) * sin (phik (j) - zeta (j));

      for (int k = 0; k < nres; k++)
	{
	  Cosk (j) += EEh (j, k) * (cos (xih (j, k)) * Xk (k) - sin (xih (j, k)) * Yk (k));
	  Sink (j) += EEh (j, k) * (cos (xih (j, k)) * Yk (k) + sin (xih (j, k)) * Xk (k));
	  sink (j) += EEh (j, k) * Psik (k) * sin (phik (j) - phik (k) - xih (j, k));
	}
    }

  Array<double,1> XkRHS      (nres);
  Array<double,1> YkRHS      (nres);
  Array<double,2> alphakpRHS (nres, NFLOW);
  Array<double,2> betakpRHS  (nres, NFLOW);

  for (int j = 0; j < nres; j++)
    {
      double omega = GetActualFrequency (j);
      double Wk    = (0.8227/2.) * 4. * fack (j) * sqrt (fabs (Psik (j))) /a (j) /rk (j);
      double dk    = delk (j) /(R_0 * a (j)) /rk (j);

      if (LIN)
	{
	  XkRHS (j) = - omega * Yk (j) + Cosk (j) /Sk (j) /dk;
	  YkRHS (j) = + omega * Xk (j) + Sink (j) /Sk (j) /dk;
	}
      else
	{
	  XkRHS (j) = - omega * Yk (j) + Cosk (j) /Sk (j) /(Wk + dk);
	  YkRHS (j) = + omega * Xk (j) + Sink (j) /Sk (j) /(Wk + dk);
	}
      
      for (int i = 0; i < NFLOW; i++)
	{
	  alphakpRHS (j, i) = (torp (j, i) * Psik (j) * sink (j)
			       - (j1p (i)*j1p (i) /taumk (j) + 1. /tautk (j) + 1. /tauxk (j)) * alphakp (j, i))
	    /(1. + 2.*qhatk (j)*qhatk (j));

	  betakpRHS (j, i) = tort (j, i) * Psik (j) * sink (j)
	    - (j0p (i)*j0p (i) /taumk (j)  + 1. /tauxk (j)) * betakp (j, i);
	}
    }

  PackRhs (XkRHS, YkRHS, alphakpRHS, betakpRHS, dydt);
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
void Phase::RK4Adaptive (double& x, Array<double,1>& y, double& h, 
			 double& t_err, double acc, double S, int& rept,
			 int maxrept, double h_min, double h_max, int flag, 
			 int diag, FILE* file)
{
  int neqns = y.extent (0);
  Array<double,1> y0 (neqns), y1 (neqns);
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
          err   = fabs (y (i) - y1 (i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs ((y (i) - y1 (i)) /y (i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = fabs ((y (i) - y1 (i)) /y (i));
          err2  = fabs (y (i) - y1 (i));
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
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
  h = (fabs (h) > h_max) ? h_max * h /fabs (h) : h;

  // Abort if step-length falls below h_min
  if (fabs (h) < h_min)
    { 
      //printf ("Phase::RK4Adpative: Warning - |h| < hmin at x = %11.4e\n", x);
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
void Phase::RK4Fixed (double& x, Array<double,1>& y, double h)
{
  int neqns = y.extent (0);
  Array<double,1> dydx (neqns), k1 (neqns), k2 (neqns), k3 (neqns);
  Array<double,1> k4   (neqns), f  (neqns);

  // Zeroth intermediate step 
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1 (i) = h * dydx (i);
      f  (i) = y (i) + k1 (i) /2.;
    }

  // First intermediate step 
  Rhs (x + h /2., f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2 (i) = h * dydx (i);
      f  (i) = y (i) + k2 (i) /2.;
    }

  // Second intermediate step 
  Rhs (x + h /2., f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3 (i) = h * dydx (i);
      f  (i) = y (i) + k3 (i);
    }

  // Third intermediate step 
  Rhs (x + h, f, dydx);
  for (int i = 0; i < neqns; i++)
    k4 (i) = h * dydx (i);

  // Actual step 
  for (int i = 0; i < neqns; i++)
    y (i) += k1 (i) /6. + k2 (i) /3. + k3 (i) /3. + k4 (i) /6.;
  x += h;
}

// #################################
// Function to open file for writing
// #################################
FILE* Phase::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("PHASE::OpenFilew: Error opening data-file %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to open file for reading
// #################################
FILE* Phase::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("PHASE::OpenFilew: Error opening data-file %s\n", filename);
      exit (1);
    }
  return file;
}

// ###################################
// Function to open file for appending
// ###################################
FILE* Phase::OpenFilea (char* filename)
{
  FILE* file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("PHASE::OpenFilea: Error opening data-file %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to call operating system
// #################################
void Phase::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("PHASE: Operating system call error executing %s\n", command);
      exit (1);
    }
}
