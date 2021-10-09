// Phase.cpp

// #####################
// PROGRAM ORGANIZATION:
// #####################

//        Phase:: Phase                     ()
// void   Phase:: Solve                     ()
// void   Phase:: Read_Data                 ()
// void   Phase:: Scan_Shift                ()
// void   Phase:: Calc_Velocity             ()
// void   Phase:: Initialize                ()
// void   Phase:: Save                      ()
// void   Phase:: IslandDynamics            ()
// void   Phase:: CalcRMP                   (double t)
// void   Phase:: CalcCoil                  (double t, double& IU, double& IM, double& IL, double& PU, double& PM, double& PL)
// int    Phase:: Findk                     ()
// double Phase:: FindMax                   (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML)
// double Phase:: FindMin                   (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML)
// void   Phase:: ThreeCoil                 (double XUM, double XUL, double XML, double gammaUM, double gammaUL, double gammaML, double Delta, double& fun, double& deriv, double& dderiv)
// void   Phase:: CalcChiZeta               (double t)
// void   Phase:: Pack                      (Array<double,1> y)
// void   Phase:: Unpack                    (Array<double,1> y)
// void   Phase:: PackRhs                   (Array<double,1> XkRHS, Array<double,1> YkRHS,
//		                             Array<double,2> alphakpRHS, Array<double,2> betakpRHS, Array<double,1> dydt)
// double Phase:: GetNaturalFrequency       (int j)
// double Phase:: GetActualFrequency        (int j)
// double Phase:: GetDeltaOmega             (int j)
// double Phase:: GetDeltaOmegaTheta        (int j)
// double Phase:: GetDeltaOmegaPhi          (int j)
// double Phase:: GetDeltaVPhi              (int j)
// double Phase:: GetDeltaVParallel         (int j)
// double Phase:: GetDeltaVEB               (int j)
// double Phase:: GetDeltaEr                (int j)
// void   Phase:: GetElectromagneticTorques (Array<double,1 >& y, double* T_Rmp, double* T_Wall, double* T_Tear)
// void   Phase:: Rhs                       (double t, Array<double,1>& y, Array<double,1>& dydt)
// void   Phase:: RK4Adaptive               (double& x, Array<double,1>& y, double& h, 
//			                     double& t_err, double acc, double S, int& rept,
//			                     int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
// void   Phase:: RK4Fixed                  (double& x, Array<double,1>& y, double h)
// FILE*  Phase:: OpenFilew                 (char* filename)
// FILE*  Phase:: OpenFiler                 (char* filename)
// FILE*  Phase:: OpenFilea                 (char* filename)
// void   Phase:: CallSystem                (char* command)

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
void Phase::Solve ()
{
  // Read input data
  Read_Data ();
  fflush (stdout);

  // Scan RMP phase shift
  Scan_Shift ();
  fflush (stdout);

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
  gsl_vector_complex_free (DeltaU);  gsl_vector_complex_free (DeltaL);
  gsl_vector_complex_free (DeltaM);  gsl_vector_complex_free (ChiU);
  gsl_vector_complex_free (ChiL);    gsl_vector_complex_free (ChiM);
}

// ###########################
// Function to read input data
// ###########################
void Phase::Read_Data ()
{
  // Output version information
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  
  // ......................................
  // Set default values of input parameters
  // ......................................
  NFLOW   = 200;

  STAGE5  = 1;

  INTF    = 0;
  INTN    = 0;
  INTU    = 0;
  NATS    = 0;
  OLD     = 0;

  LIN     = 0;
  FREQ    = 0;
  FFAC    = 0.;

  CXD     = 1;
  BOOT    = 1;
  CURV    = 1;
  POLZ    = 1;

  MID     = 2;
  COPT    = 0;
  CORE    = 1.;

  DT      = 1.e-5;
  TSTART  = 0.;
  TEND    = 1.e6;

  SCALE   = 2.;
  PMAX    = 4.;
  NPHA    = 361;
  CHIR    = 1.;
  HIGH    = 1;
  RATS    = 1;

  TAUW    = 1.;
  TOFF    = 0.;

  TYPE    = 1;
  SSTART  = 10.;
  SEND    = 15.;
  WAMOD   = 0.;
  WPMOD   = 0.;
  SAMP    = 1.;
  SPHA    = 0.;
  BACK    = 0.1;
  RPERIOD = 1000.;
  RSTART  = 0.;
  REND    = 4.;
  RPHA    = 0.5;
  
  // Read input data from namelists (Inputs/Phase.nml, Inputs/Waveform.nml)
  printf ("........................................................................................\n");
  printf ("Input parameters (from Inputs/Phase.nml, Inputs/Waveform.nml, and command line options):\n");
  printf ("........................................................................................\n");

  TCTRL = new double[MAXCONTROLPOINTNUMBER];
  ICTRL = new double[MAXCONTROLPOINTNUMBER];
  PCTRL = new double[MAXCONTROLPOINTNUMBER];

  NameListRead (&NFLOW, &STAGE5, &INTF, &INTN, &INTU, &NATS, &OLD, &FREQ, &LIN, &MID, &COPT,
		&DT, &TSTART, &TEND, &TOFF, &SCALE, &PMAX, &CHIR, &HIGH, &RATS, &CORE, &FFAC, 
		&CXD, &BOOT, &CURV, &POLZ, &TAUW, &TYPE, &NCTRL, TCTRL, ICTRL, PCTRL,
		&SSTART, &SEND, &WAMOD, &WPMOD, &SAMP, &SPHA, &BACK, &RPERIOD, &RSTART, &REND, &RPHA);

  TT.resize (NCTRL);
  
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
  if (FREQ < 0 || FREQ > 3)
    {
      printf ("PHASE:: Invalid FREQ value\n");
      exit (1);
    }
  if (MID < 1 || MID > 3)
    {
      printf ("PHASE:: Invalid MID value\n");
      exit (1);
    }
  if (TAUW <= 0.)
    {
      printf ("PHASE:: Invalid TAUW value\n");
      exit (1);
    }
  if (TOFF < 0.)
    {
      printf ("PHASE:: Invalid TOFF value\n");
      exit (1);
    }
  if (SSTART > SEND)
    {
      printf ("PHASE:: SSTART must be less than SEND\n");
      exit (1);
    }
  if (TYPE < 0 || TYPE > 1)
    {
      printf ("PHASE:: Invalid TYPE value\n");
      exit (1);
    }
  if (RPERIOD <= 0.)
    {
      printf ("PHASE:: Invalid RPERIOD value\n");
      exit (1);
    }
  
  for (int i = 0; i < NCTRL; i++)
    printf ("T = %10.3e  IRMP = %10.3e  PRMP/pi = %10.3e\n", TCTRL[i], ICTRL[i], PCTRL[i]/M_PI);

  // .............................
  // Output calculation parameters
  // .............................
  printf ("NFLOW = %4d STAGE5 = %2d INTF = %2d INTN = %2d INTU = %2d OLD = %2d LIN = %2d MID = %2d COPT = %2d CORE = %10.3e HIGH = %2d RATS = %2d\n",
	  NFLOW, STAGE5, INTF, INTN, INTU, OLD, LIN, MID, COPT, CORE, HIGH, RATS);
  printf ("FREQ = %2d FFAC = %10.3e CXD = %2d BOOT = %2d CURV = %2d POLZ = %2d TAUW = %10.3e DT = %10.3e TSTART = %10.3e TEND = %10.3e TOFF = %10.3e SCALE = %10.3e PMAX = %10.3e CHIR = %10.3e\n",
	  FREQ, FFAC, CXD, BOOT, CURV, POLZ, TAUW, DT, TSTART, TEND, TOFF, SCALE, PMAX, CHIR);
  printf ("TYPE = %2d NCTRL = %4d SSTART = %10.3e SEND = %10.3e WAMOD = %10.3e WPMOD = %10.3e SAMP = %10.3e SPHA = %10.3e BACK = %10.3e RPERIOD = %10.3e RSTART = %10.3e REND = %10.3e RPHA = %10.3e\n",
	  TYPE, NCTRL, SSTART, SEND, WAMOD, WPMOD, SAMP, SPHA, BACK, RPERIOD, RSTART, REND, RPHA);
  
  FILE* namelist = OpenFilew ((char*) "Outputs/InputParameters.txt");
  fprintf (namelist, "Git Hash     = "); fprintf (namelist, GIT_HASH);     fprintf (namelist, "\n");
  fprintf (namelist, "Compile time = "); fprintf (namelist, COMPILE_TIME); fprintf (namelist, "\n");
  fprintf (namelist, "Git Branch   = "); fprintf (namelist, GIT_BRANCH);   fprintf (namelist, "\n\n");
  fprintf (namelist, "Input parameters (from Inputs/Phase.nml and command line options):\n");
  fprintf (namelist, "NFLOW = %4d STAGE5 = %2d INTF = %2d INTN = %2d INTU = %2d NATS = %2d OLD = %2d LIN = %2d MID = %2d COPT = %2d CORE = %10.3e HIGH = %2d RATS = %2d \n",
	   NFLOW, STAGE5, INTF, INTN, INTU, NATS, OLD, LIN, MID, COPT, CORE, HIGH, RATS);
  fprintf (namelist, "FREQ = %2d FFAC = %10.3e CXD = %2d BOOT = %2d CURV = %2d POLZ = %2d TAUW = %10.3e DT = %10.3e TSTART = %10.3e TEND = %10.3e TOFF = %10.3e SCALE = %10.3e PMAX = %10.3e CHIR = %10.3e\n",
	   FREQ, FFAC, CXD, BOOT, CURV, POLZ, TAUW, DT, TSTART, TEND, TOFF, SCALE, PMAX, CHIR);
  fprintf (namelist, "TYPE = %2d NCTRL = %4d SSTART = %10.3e SEND = %10.3e WAMOD = %10.3e WPMOD = %10.3e SAMP = %10.3e SPHA = %10.3e BACK = %10.3e RPERIOD = %10.3e RSTART = %10.3e REND = %10.3e RPHA = %10.3e\n",
	  TYPE, NCTRL, SSTART, SEND, WAMOD, WPMOD, SAMP, SPHA, BACK, RPERIOD, RSTART, REND, RPHA);
  fclose (namelist);

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
   
  // ...........................
  // Read data from program FLUX
  // ...........................
  printf ("...............................\n");
  printf ("Reading data from program FLUX:\n");
  printf ("...............................\n");
  
  int ini; double inr; int NPSI;
  double Freal, Fimag;
  double Ereal, Eimag; 
 
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
	  printf ("PHASE: Error reading fFile (2)\n");
	  exit (1);
	}
    }

  A1.resize     (nres);
  qhatk.resize  (nres);
  C1k.resize    (nres);
  C2k.resize    (nres);
  Poem1.resize  (nres);
  Poem2.resize  (nres);
  Poem3.resize  (nres);
  Deltaw.resize (nres);
  Sigmaw.resize (nres);
  for (int j = 0; j < nres; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &ini, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &A1(j), &inr, &qhatk(j), &C1k(j), &C2k(j), &inr,
		  &Poem1(j), &Poem2(j), &Poem3(j), &Deltaw(j), &Sigmaw(j)) != 23)
	{
	  printf ("PHASE: Error reading fFile (3)\n");
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
 
  fclose (file);

  printf ("R_0 = %10.3e  B_0 = %10.3e  nres = %3d\n", R_0, B_0, nres);

  printf ("E-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%10.3e,%10.3e) ", GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
      printf ("\n");
    }
 
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

  printf ("tau_A = %10.3e  P0 = %10.3e\n", tau_A, P0);

  // Normalize times
  for (int i = 0; i < NCTRL; i++)
    TT (i) = TCTRL[i] *1.e-3/tau_A;
  dTT     = DT        *1.e-3/tau_A;
  Tstart  = TSTART    *1.e-3/tau_A;
  Tend    = TEND      *1.e-3/tau_A;
  Toff    = TOFF      *1.e-3/tau_A;
  tauw    = TAUW      *1.e-3/tau_A;
  SSTART  = SSTART    *1.e-3/tau_A;
  SEND    = SEND      *1.e-3/tau_A;
  RPERIOD = RPERIOD   *1.e-3/tau_A;

  // Normalize frequencies
  WAMOD = WAMOD * tau_A/1.e-3;
  WPMOD = WPMOD * tau_A/1.e-3;

  mk.resize      (nres); ntor.resize     (nres); rk.resize       (nres); qk.resize       (nres); rhok.resize   (nres);
  a.resize       (nres); Sk.resize       (nres); taumk.resize    (nres); tautk.resize    (nres); tauxk.resize  (nres);
  fack.resize    (nres); delk.resize     (nres); wkl.resize      (nres); wke.resize      (nres); wkn.resize    (nres);
  dnedrk.resize  (nres); dTedrk.resize   (nres); Wcrnek.resize   (nres); WcrTek.resize   (nres); WcrTik.resize (nres);
  akk.resize     (nres); gk.resize       (nres); dPsiNdr.resize  (nres); PsiN.resize     (nres); nek.resize    (nres);
  nik.resize     (nres); Tek.resize      (nres); Tik.resize      (nres); dnidrk.resize   (nres); dTidrk.resize (nres);
  Factor1.resize (nres); Factor2.resize  (nres); Factor3.resize  (nres); Factor4.resize  (nres); 
  Factor5.resize (nres); Factor6.resize  (nres); Factor7.resize  (nres); Factor8.resize  (nres);
  Factor9.resize (nres); Factor10.resize (nres); Factor11.resize (nres); Factor12.resize (nres);
  alphabe.resize (nres); alphabi.resize  (nres); alphac.resize   (nres); alphap.resize   (nres);
  rhothe.resize  (nres); rhothi.resize   (nres); etae.resize     (nres); etai.resize     (nres); chipk.resize  (nres);

  for (int j = 0; j < nres; j++)
    if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&mk      (j), &ntor     (j), &rk       (j), &qk       (j), &rhok   (j),
		&a       (j), &Sk       (j), &taumk    (j), &tautk    (j), &fack   (j),
		&delk    (j), &wkl      (j), &wke      (j), &wkn      (j),
		&dnedrk  (j), &dTedrk   (j), &Wcrnek   (j), &WcrTek   (j), &WcrTik (j),
		&akk     (j), &gk       (j), &dPsiNdr  (j), &PsiN     (j), &nek    (j),
		&nik     (j), &Tek      (j), &Tik      (j), &dnidrk   (j), &dTidrk (j),
		&Factor1 (j), &Factor2  (j), &Factor3  (j), &Factor4  (j),
		&Factor5 (j), &Factor6  (j), &Factor7  (j), &Factor8  (j),
		&Factor9 (j), &Factor10 (j), &Factor11 (j), &Factor12 (j), &tauxk (j),
		&alphabe (j), &alphabi  (j), &alphac   (j), &alphap   (j), &rhothe (j),
		&rhothi  (j), &etae     (j), &etai     (j), chipk     (j)) != 51)
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
    printf ("m = %3d h_r = %9.3e q = %9.3e g = %9.3e akk = %9.3e h_rho = %9.3e h_a = %9.3e S = %9.3e h_tauM = %9.3e h_tauth = %9.3e h_del = %9.3e A1 = %9.3e q_hat = %9.3e abe = %10.3e abi = %10.3e ac = %10.3e ap = %10.3e rthe = %9.3e rthi = %9.3e\n",
	    mk (j), rk (j), qk (j), gk (j), akk (j), rhok (j), a (j), Sk (j), taumk (j), tautk (j), delk (j) /(rk (j) * a (j) * R_0), A1 (j), qhatk (j),
	    alphabe (j), alphabi (j), alphac (j), alphap (j), rhothe (j), rhothi (j));

  // Set Deltak+/- values
  Deltakp.resize (nres); Deltakm.resize (nres);

  Deltakm (0) = CHIR * PsiN (0);
  for (int j = 1; j < nres; j++)
    Deltakm (j) = CHIR * (PsiN (j) - PsiN (j-1));

  for (int j = 0; j < nres-1; j++)
    Deltakp (j) = CHIR * (PsiN (j+1) - PsiN (j));
  Deltakp (nres-1) = CHIR * (PSILIM - PsiN (nres-1));

  // Implement CXD, BOOT, CURV, and POLZ flags
  for (int j = 0; j < nres; j++)
    {
      if (CXD == 0)
	tauxk (j) = 1.e15;
      if (BOOT == 0)
	{
	  alphabe (j) = 0.;
	  alphabi (j) = 0.;
	}
      if (CURV == 0)
	alphac (j) = 0.;
      if (POLZ == 0)
	alphap (j) = 0.;
    }
  
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
  
  char    line[MAXULFILELINELENGTH];
  char    line1[MAXULFILELINELENGTH];
  char    line0[MAXULFILELINELENGTH];
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
  printf ("SCALEFACTOR = %10.3e\n", SCALEFACTOR);

  if (MID == 3)
    {
      printf ("Middle coil:\n");
      file = OpenFiler ((char*) "Inputs/mFile");

      fgets (line0, MAXULFILELINELENGTH, file);
      for (int i = 0; i < 4; i++)
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

	    if (strstr (line0, "GPEC") != NULL)
	      {
		DRE[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
		DIM[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
	      }
	    else
	      {
		DRE[i] /= 2.*M_PI * qk(i);
		DIM[i] /= 2.*M_PI * qk(i);
	      }
	    
	    CRE[i] =   DIM[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    CIM[i] = - DRE[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    
	    gsl_vector_complex_set (DeltaM, i, gsl_complex_rect (DRE[i], DIM[i]));
	    gsl_vector_complex_set (ChiM,   i, gsl_complex_rect (CRE[i], CIM[i]));
	    
	    double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiM, i));
	    double WUNRE = 4. * sqrt (A1 (i) * Psi);
	    double WFULL = sqrt (FFh (i, i) * EEh (i, i)) * WUNRE;

	    printf ("GPEC: q = %10.3e  Psi = %10.3e  PsiN = %10.3e  Delta = (%10.3e, %10.3e)  Chi = (%10.3e, %10.3e)  W_UNRE = %10.3e  W_UNRE/W_GPEC = %10.3e  W_FULL/W_GPEC = %10.3e\n",
		    QIN[i], PSI[i], PsiN(i), DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
	  }
      fclose (file);

      if (fabs (QIN[0] - qk(0)) > 1.e-3)
	{
	  printf ("PHASE:: Error - minimum resonant q values do not match in nFile and mFile\n");
	  exit (1);
	}
    }

  if (MID >= 2)
    {
      printf ("Upper coil:\n");
      file = OpenFiler ((char*) "Inputs/uFile");

      fgets (line0, MAXULFILELINELENGTH, file);
      for (int i = 0; i < 4; i++)
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

	    if (strstr (line0, "GPEC") != NULL)
	      {
		DRE[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
		DIM[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
	      }
	    else
	      {
		DRE[i] /= 2.*M_PI * qk(i);
		DIM[i] /= 2.*M_PI * qk(i);
	      }
	    
	    CRE[i] =   DIM[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk (i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    CIM[i] = - DRE[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk (i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    
	    gsl_vector_complex_set (DeltaU, i, gsl_complex_rect (DRE[i], DIM[i]));
	    gsl_vector_complex_set (ChiU,   i, gsl_complex_rect (CRE[i], CIM[i]));
	    
	    double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiU, i));
	    double WUNRE = 4. * sqrt (A1 (i) * Psi);
	    double WFULL = sqrt (FFh (i, i) * EEh (i, i)) * WUNRE;
	    
	    printf ("GPEC: q = %10.3e  Psi = %10.3e  PsiN = %10.3e  Delta = (%10.3e, %10.3e)  Chi = (%10.3e, %10.3e)  W_UNRE = %10.3e  W_UNRE/W_GPEC = %10.3e  W_FULL/W_GPEC = %10.3e\n",
		    QIN[i], PSI[i], PsiN(i), DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
	  }
      fclose (file);

      if (fabs (QIN[0] - qk(0)) > 1.e-3)
	{
	  printf ("PHASE:: Error - minimum resonant q values do not match in nFile and uFile\n");
	  exit (1);
	}
    }

  if (MID >= 1)
    {
      printf ("Lower coil:\n");
      file = OpenFiler ((char*) "Inputs/lFile");

      fgets (line0, MAXULFILELINELENGTH, file);
      for (int i = 0; i < 4; i++)
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

	    if (strstr (line0, "GPEC") != NULL)
	      {
		DRE[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
		DIM[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
	      }
	    else
	      {
		DRE[i] /= 2.*M_PI * qk(i);
		DIM[i] /= 2.*M_PI * qk(i);
	      }
	    
	    CRE[i] =   DIM[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    CIM[i] = - DRE[i] * (rk (i) * a (i)) * (rk (i) * a (i)) * gk(i)
	      /double (mk (i)) /(akk (i) + rk (i) * rk (i) * a (i) * a (i) /qk (i) /qk (i)) /EEh (i, i);
	    
	    gsl_vector_complex_set (DeltaL, i, gsl_complex_rect (DRE[i], DIM[i]));
	    gsl_vector_complex_set (ChiL,   i, gsl_complex_rect (CRE[i], CIM[i]));
	    
	    double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiL, i));
	    double WUNRE = 4. * sqrt (A1 (i) * Psi);
	    double WFULL = sqrt (FFh (i, i) * EEh (i, i)) * WUNRE;

	    printf ("GPEC: q = %10.3e  Psi = %10.3e  PsiN = %10.3e  Delta = (%10.3e, %10.3e)  Chi = (%10.3e, %10.3e)  W_UNRE = %10.3e  W_UNRE/W_GPEC = %10.3e  W_FULL/W_GPEC = %10.3e\n",
		    QIN[i], PSI[i], PsiN(i), DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
	  }
      fclose (file);

      if (fabs (QIN[0] - qk(0)) > 1.e-3)
	{
	  printf ("PHASE:: Error - minimum resonant q values do not match in nFile and lFile\n");
	  exit (1);
	}
    } 
 
  delete[] QIN; delete[] DRE; delete[] DIM; delete[] CRE; delete[] CIM; delete[] WWW; delete[] PSI;
}

// ########################################################################
// Function to calculate vacuum flux versus relative RMP coil current phase
// ########################################################################
void Phase::Scan_Shift ()
{
  // Perform scan
  phase.resize (NPHA);
  wvac.resize  (nres, NPHA);
  double one = 1.;

  for (int i = 0; i < NPHA; i++)
    {
      double      pha   = - PMAX*M_PI/2. + double (i) * PMAX*M_PI /double (NPHA-1);
      phase (i)         = pha;
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

	  double chi  =   gsl_complex_abs (h);
	  double zeta = - gsl_complex_arg (h);
	  double wv   = 4. * R_0 * fack (j) * sqrt (chi);

	  wvac (j, i) = GetVacuumIslandWidth (j, chi);
	}
    }

  // Write Stage 4 netcdf file
  WriteStage4Netcdfc ();
}

// ####################################
// Function to write Stage4 NETCDF file
// ####################################
void Phase::WriteStage4Netcdfc ()
{
  // Convert data from blitz++ array to c array
  int* mk_x = new int[nres];
  for (int i = 0; i < nres; i++)
    mk_x[i] = mk(i);
  double* phase_x = new double[NPHA];
  for (int i = 0; i < NPHA; i++)
    phase_x[i] = phase(i) /M_PI;

  // Open file
  int err = 0, dataFile;
  err = nc_create ("Outputs/Stage4.nc", NC_CLOBBER, &dataFile);
  
  if (err != 0)
    {
      printf ("PHASE::WriteStage2Netcdfc: Error opening Outputs/Stage4.nc\n");
      exit (1);
    }
  
  // m_k
  int nres_d, mk_y;
  err += nc_def_dim (dataFile, "i_res", nres, &nres_d);
  err += nc_def_var (dataFile, "m_pol", NC_INT, 1, &nres_d, &mk_y);

  // phase
  int P_d, phase_y;
  err += nc_def_dim (dataFile, "i_phase", NPHA, &P_d);
  err += nc_def_var (dataFile, "phase",   NC_DOUBLE, 1, &P_d, &phase_y);
 
  // w_vac
  double* DATA = new double[nres*NPHA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NPHA; j++)
      DATA[j + i*NPHA] = wvac (i, j);

  int W_d[2], wvac_y; 
  W_d[0] = nres_d;
  W_d[1] = P_d;
  err += nc_def_var (dataFile, "W_vacuum", NC_DOUBLE, 2, W_d, &wvac_y);

  err += nc_enddef (dataFile);

  if (err != 0)
    {
      printf ("PHASE::WriteStage2Netcdfc: Error defining variables in Outputs/Stage2.nc\n");
      exit (1);
    }

  // Write data
  err += nc_put_var_int    (dataFile, mk_y,    mk_x);
  err += nc_put_var_double (dataFile, phase_y, phase_x);
  err += nc_put_var_double (dataFile, wvac_y,  DATA);

  if (err != 0)
    {
      printf ("PHASE::WriteStage2Netcdfc: Error writing Outputs/Stage4.nc\n");
      exit (1);
    }
  
  // Close file
  err += nc_close (dataFile);

  if (err != 0)
    {
      printf ("PHASE::WriteStage2Netcdfc: Error closing Outputs/Stage4.nc\n");
      exit (1);
    }

  // Clean up
  delete[] mk_x; delete[] DATA; delete[] phase_x;
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
      Psik  (j) = 0.;
      phik  (j) = 0.;
      Xk    (j) = 0.;
      Yk    (j) = 0.;
      Psiwk (j) = 0.;
      phiwk (j) = 0.;
      Uk    (j) = 0.;
      Vk    (j) = 0.;
 
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
      Array<double,1> _Psiwk   (_nres);
      Array<double,1> _phiwk   (_nres);
      Array<double,2> _alphakp (_nres, _NFLOW);
      Array<double,2> _betakp  (_nres, _NFLOW);
      Array<int,1>    _lock    (_nres);
      Array<double,1> _ww      (_nres);

      int in;
      for (int j = 0; j < _nres; j++)
	if (fscanf (file, "%d %lf %lf %lf %lf %d %lf\n", &in, &_Psik (j), &_phik (j), &_Psiwk (j), &_phiwk (j), &_lock (j), &_ww (j)) != 7)
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
	  Psik  (j) = _Psik (j);
	  phik  (j) = _phik (j);
	  Xk    (j) = Psik (j) * cos (phik (j));
	  Yk    (j) = Psik (j) * sin (phik (j));
	  Psiwk (j) = _Psiwk (j);
	  phiwk (j) = _phiwk (j);
	  Uk    (j) = Psiwk (j) * cos (phiwk (j));
	  Vk    (j) = Psiwk (j) * sin (phiwk (j));
 	  lock  (j) = _lock (j);
	  ww    (j) = _ww   (j);
	}

      for (int j = _nres; j < nres; j++)
	{
	  Psik  (j) = 0.;
	  phik  (j) = 0.;
	  Xk    (j) = 0.;
	  Yk    (j) = 0.;
	  Psiwk (j) = 0.;
	  phiwk (j) = 0.;
	  Uk    (j) = 0.;
	  Vk    (j) = 0.;
 	  lock  (j) = 0.;
	  ww    (j) = GetActualFrequency (j);
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
    fprintf (file, "%3d %19.6e %19.6e %19.6e %19.6e %2d %19.6e\n", j, Psik (j), phik (j), Psiwk (j), phiwk (j), lock (j), ww (j));

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
  Psiwk.resize   (nres);
  phiwk.resize   (nres);
  Uk.resize      (nres);
  Vk.resize      (nres);
  alphakp.resize (nres, NFLOW);
  betakp.resize  (nres, NFLOW);
  lock.resize    (nres);
  ww.resize      (nres);
  vp.resize      (nres);

  Initialize ();

  // .............................
  // Initialize Netcdf data arrays
  // .............................
  int*    mk_x   = new int[nres];
  double* PsiN_x = new double[nres];
  for (int i = 0; i < nres; i++)
    {
      mk_x[i]   = mk(i);
      PsiN_x[i] = PsiN(i);
    }

  int     zcount  = -1;  
  double* zTime   = new double[1];   
  double* zIrmp   = new double[1];   
  double* zPrmp   = new double[1];
  double* zIU     = new double[1];   
  double* zIM     = new double[1];
  double* zIL     = new double[1];
  double* zPU     = new double[1];   
  double* zPM     = new double[1];
  double* zPL     = new double[1];   
  double* zphi    = new double[nres];    
  double* zvph    = new double[nres];    
  double* zomega0 = new double[nres]; 
  double* zomega  = new double[nres];  
  double* zPsim   = new double[nres];   
  double* zPsip   = new double[nres];   
  double* zW      = new double[nres];      
  double* zPsivm  = new double[nres];  
  double* zPsivp  = new double[nres];  
  double* zWv     = new double[nres];     
  double* zdP     = new double[nres];     
  double* zdP0    = new double[nres];
  double* zTRmp   = new double[nres];
  double* zTWall  = new double[nres];
  double* zTTear  = new double[nres];
  double* ztlock  = new double[nres];
  double* zIlock  = new double[nres];
  double* znelock = new double[nres];
  double* zTelock = new double[nres];
  double* zCplock = new double[nres];
  double* zB0lock = new double[nres];
  double* zR0lock = new double[nres];
  double* zwllock = new double[nres];
  double* zwelock = new double[nres];
  double* zwnlock = new double[nres];
  double* zbrlock = new double[nres];
  double* zWlock  = new double[nres];

  for (int i = 0; i < nres; i++)
    {
      ztlock[i]  = -1000.;
      zIlock[i]  = 0.;
      znelock[i] = 0.;
      zTelock[i] = 0.;
      zCplock[i] = 0.;
      zB0lock[i] = 0.;
      zR0lock[i] = 0.;
      zwllock[i] = 0.;
      zwelock[i] = 0.;
      zwnlock[i] = 0.;
      zbrlock[i] = 0.;
      zWlock[i]  = 0.;
    }

  size_t* tstart = new size_t[1];
  size_t* tcount = new size_t[1];
  size_t* dstart = new size_t[2];
  size_t* dcount = new size_t[2];

  tcount[0] = 1;
  dcount[0] = 1;
  dstart[1] = 0;
  dcount[1] = nres;

  // ......................
  // Initialize Netcdf file
  // ......................

  int err, dataFile;
  err = nc_create ("Outputs/Stage5.nc", NC_CLOBBER, &dataFile);

  if (err != 0)
    {
      printf ("PHASE::IslandDynamics: Error opening Outputs/Stage5.nc\n");
      exit (1);
    }

  int nres_d, mk_y;
  err += nc_def_dim (dataFile, "i_res", nres, &nres_d);
  err += nc_def_var (dataFile, "m_pol", NC_INT, 1, &nres_d, &mk_y);

  int PsiN_y;
  err += nc_def_var (dataFile, "PsiN_res", NC_DOUBLE, 1, &nres_d, &PsiN_y);

  int time_d, dim[2];
  err += nc_def_dim (dataFile, "time",  NC_UNLIMITED, &time_d);
  dim[0] = time_d;
  dim[1] = nres_d;

  int yTime;
  err += nc_def_var (dataFile, "time",     NC_DOUBLE, 1, &time_d, &yTime);

  int yIrmp;
  err += nc_def_var (dataFile, "Irmp",     NC_DOUBLE, 1, &time_d, &yIrmp);

  int yPrmp;
  err += nc_def_var (dataFile, "Prmp",     NC_DOUBLE, 1, &time_d, &yPrmp);

  int yIU;
  err += nc_def_var (dataFile, "IU",       NC_DOUBLE, 1, &time_d, &yIU);

  int yIM;
  err += nc_def_var (dataFile, "IM",       NC_DOUBLE, 1, &time_d, &yIM);

  int yIL;
  err += nc_def_var (dataFile, "IL",       NC_DOUBLE, 1, &time_d, &yIL);

  int yPU;
  err += nc_def_var (dataFile, "PU",       NC_DOUBLE, 1, &time_d, &yPU);

  int yPM;
  err += nc_def_var (dataFile, "PM",       NC_DOUBLE, 1, &time_d, &yPM);

  int yPL;
  err += nc_def_var (dataFile, "PL",       NC_DOUBLE, 1, &time_d, &yPL);

  int yphi;
  err += nc_def_var (dataFile, "phi",      NC_DOUBLE, 2, dim,     &yphi);

  int yvph;
  err += nc_def_var (dataFile, "phi_dot",  NC_DOUBLE, 2, dim,     &yvph);

  int yomega0;
  err += nc_def_var (dataFile, "omega0",   NC_DOUBLE, 2, dim,     &yomega0);

  int yomega;
  err += nc_def_var (dataFile, "omega",    NC_DOUBLE, 2, dim,     &yomega);

  int yW;
  err += nc_def_var (dataFile, "W",        NC_DOUBLE, 2, dim,     &yW);

  int yWv;
  err += nc_def_var (dataFile, "W_vac",    NC_DOUBLE, 2, dim,     &yWv);

  int yPsim;
  err += nc_def_var (dataFile, "Psi_-",    NC_DOUBLE, 2, dim,     &yPsim);

  int yPsip;
  err += nc_def_var (dataFile, "Psi_+",    NC_DOUBLE, 2, dim,     &yPsip);
  
  int yPsivm;
  err += nc_def_var (dataFile, "Psi_v-",   NC_DOUBLE, 2, dim,     &yPsivm);

  int yPsivp;
  err += nc_def_var (dataFile, "Psi_v+",   NC_DOUBLE, 2, dim,     &yPsivp);

  int ydP;
  err += nc_def_var (dataFile, "DeltaP",   NC_DOUBLE, 2, dim,     &ydP);

  int ydP0;
  err += nc_def_var (dataFile, "DeltaP0",  NC_DOUBLE, 2, dim,     &ydP0);

  int yTRmp;
  err += nc_def_var (dataFile, "T_Rmp",    NC_DOUBLE, 2, dim,     &yTRmp);

  int yTWall;
  err += nc_def_var (dataFile, "T_Wall",   NC_DOUBLE, 2, dim,     &yTWall);

  int yTTear;
  err += nc_def_var (dataFile, "T_Tear",   NC_DOUBLE, 2, dim,     &yTTear);

  int ytlock;
  err += nc_def_var (dataFile, "t_lock",   NC_DOUBLE, 1, &nres_d, &ytlock);

  int yIlock;
  err += nc_def_var (dataFile, "I_lock",   NC_DOUBLE, 1, &nres_d, &yIlock);

  int ynelock;
  err += nc_def_var (dataFile, "ne_lock",  NC_DOUBLE, 1, &nres_d, &ynelock);

  int yTelock;
  err += nc_def_var (dataFile, "Te_lock",  NC_DOUBLE, 1, &nres_d, &yTelock);

  int yCplock;
  err += nc_def_var (dataFile, "Cp_lock",  NC_DOUBLE, 1, &nres_d, &yCplock);

  int yB0lock;
  err += nc_def_var (dataFile, "B0_lock",  NC_DOUBLE, 1, &nres_d, &yB0lock);

  int yR0lock;
  err += nc_def_var (dataFile, "R0_lock",  NC_DOUBLE, 1, &nres_d, &yR0lock);

  int ywllock;
  err += nc_def_var (dataFile, "wl_lock",  NC_DOUBLE, 1, &nres_d, &ywllock);

  int ywelock;
  err += nc_def_var (dataFile, "we_lock",  NC_DOUBLE, 1, &nres_d, &ywelock);

  int ywnlock;
  err += nc_def_var (dataFile, "wn_lock",  NC_DOUBLE, 1, &nres_d, &ywnlock);

  int ybrlock;
  err += nc_def_var (dataFile, "br_lock",  NC_DOUBLE, 1, &nres_d, &ybrlock);

  int yWlock;
  err += nc_def_var (dataFile, "W_lock",   NC_DOUBLE, 1, &nres_d, &yWlock);

  err += nc_enddef (dataFile);

  if (err != 0)
    {
      printf ("PHASE::IslandDynamics: Error defining variables in Outputs/Stage5.nc\n");
      exit (1);
    }

  err += nc_put_var_int    (dataFile, mk_y,   mk_x);
  err += nc_put_var_double (dataFile, PsiN_y, PsiN_x);
    
  // Integrate equations of motion
  chi.resize  (nres);
  zeta.resize (nres);

  double          t, h, t_err;
  int             rept, step = 0, skip = 0; count = 0;
  Array<double,1> y    (2*nres*(2+NFLOW));
  Array<double,1> dydt (2*nres*(2+NFLOW));
  Array<double,1> wwo  (nres);

  if (OLD == 0)
    t = Tstart - Toff;
  else
    t = Tstart;
  h = h0;
  Pack (y);
 
  int    cnt    = 0;
  double dt     = dTT - 1.e-15; 
  fflush (stdout);
  
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

      if (OLD == 0)
	{
	  if (t > Tstart)
	    dt += h;
	  else
	    dt += 0.;
	}
      else
	dt += h;
      
      // Update time in ms
      TIME = t * tau_A * 1.e3;

      // Output Stage6 mode locking data for IslandDynamics
      for (int j = 0; j < nres; j++)
	ww (j) = GetActualFrequency (j);
      
      for (int j = 0; j < nres; j++)
	if (ww (j) * wwo (j) < 0. && lock (j) == 0)
	  {
	    CalcRMP     (t);
	    CalcChiZeta (t);

	    printf ("m = %3d locks at t = %10.3e s  irmp = %10.3e kA  prmp/pi = %10.3e\n",
		    mk (j), t*tau_A, irmp, prmp /M_PI);

	    ztlock[j]  = t *tau_A/1.e-3;
	    zIlock[j]  = irmp;
	    znelock[j] = nek(j);
	    zTelock[j] = Tek(j);
	    zCplock[j] = chipk(j);
	    zB0lock[j] = B_0;
	    zR0lock[j] = R_0;
	    zwllock[j] = wkl(j) /tau_A/1.e3;
	    zwelock[j] = wke(j) /tau_A/1.e3;
	    zwnlock[j] = wkl(j) /tau_A/1.e3;
	    zbrlock[j] = double (mk(j)) * chi (j) /rk(j) /a(j);
	    zWlock[j]  = GetIslandWidth (j);
	 
	    lock (j) = 1;
	  }

      // Output Stage5 data
      if (dt > dTT)
	{
	  zcount += 1;
	  dt = 0.;

	  zTime[0] = t*tau_A/1.e-3;

	  // Output phases of reconnected fluxes
	  for (int j = 0; j < nres; j++)
	    {
	      double phase = atan2 (sin (phik (j) - zeta (j)), cos (phik (j) - zeta (j))) /M_PI;
	      zphi[j] = phase;
	    }

	  // Output modified natural frequencies
	  for (int j = 0; j < nres; j++)
	    {
	      zomega[j] = ww (j) /tau_A/1.e3;
	    }

	  // Output unmodified natural frequencies
	  for (int j = 0; j < nres; j++)
	    {
	      zomega0[j] = GetNaturalFrequency (j) /tau_A/1.e3;
	    }

	  // Output island phase velocities
	  for (int j = 0; j < nres; j++)
	    {
	      zvph[j] = vp (j) /tau_A/1.e3;
	    }
	  
	  // Output RMP data
	  CalcRMP (t);
	  zIrmp[0] = irmp;
	  zPrmp[0] = prmp /M_PI;

	  double IU, IM, IL, PU, PM, PL;
	  CalcCoil (t, IU, IM, IL, PU, PM, PL);
	  zIU[0] = IU;
	  zIM[0] = IM;
	  zIL[0] = IL;
	  zPU[0] = PU /M_PI;
	  zPM[0] = PM /M_PI;
	  zPL[0] = PL /M_PI;

	  // Output island data
	  for (int j = 0; j < nres; j++)
	    {
	      // Calculate island width in PsiN
	      double Wpk = GetIslandWidth (j);
	      zW[j] = Wpk;
	     
	      double Xminus, Xplus;
	      GetIslandLimits (j, Psik (j), Xminus, Xplus);
	      zPsim[j] = PsiN(j) + Xminus;
	      zPsip[j] = PsiN(j) + Xplus;

	      // Calculate island width in r
	      double Wrk = R_0 * Wpk /dPsiNdr (j);

	      // Calculate vacuum island width in PsiN
	      double Wvk = GetVacuumIslandWidth (j);
	      zWv[j] = Wvk;
	      
	      GetIslandLimits (j, chi (j), Xminus, Xplus);
	      zPsivm[j] = PsiN(j) + Xminus;
	      zPsivp[j] = PsiN(j) + Xplus;

	      // Calculate density and temperature flattening widths in r
	      double deltanek = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + Wcrnek (j) * Wcrnek (j));
	      double deltaTek = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + WcrTek (j) * WcrTek (j));
	      double deltaTik = (2./M_PI) * Wpk *Wrk*Wrk /(Wrk*Wrk + WcrTik (j) * WcrTik (j));

	      // Calculate pressure reduction
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
	      zdP[j]  = deltapk /P0/Pped;
	      zdP0[j] = deltapk /P0;

	      // Calculate electromagnetic torques
	      GetElectromagneticTorques (y, zTRmp, zTWall, zTTear);

	      // Output Netcdf data
	      tstart[0] = size_t (zcount);
	      dstart[0] = size_t (zcount);

	      err += nc_put_vara_double (dataFile, yTime,   tstart, tcount, zTime);
	      err += nc_put_vara_double (dataFile, yIrmp,   tstart, tcount, zIrmp);
	      err += nc_put_vara_double (dataFile, yPrmp,   tstart, tcount, zPrmp);
	      err += nc_put_vara_double (dataFile, yIU,     tstart, tcount, zIU);
	      err += nc_put_vara_double (dataFile, yIM,     tstart, tcount, zIM);
	      err += nc_put_vara_double (dataFile, yIL,     tstart, tcount, zIL);
	      err += nc_put_vara_double (dataFile, yPU,     tstart, tcount, zPU);
	      err += nc_put_vara_double (dataFile, yPM,     tstart, tcount, zPM);
	      err += nc_put_vara_double (dataFile, yPL,     tstart, tcount, zPL);
	      err += nc_put_vara_double (dataFile, yphi,    dstart, dcount, zphi);
	      err += nc_put_vara_double (dataFile, yvph,    dstart, dcount, zvph);
	      err += nc_put_vara_double (dataFile, yomega0, dstart, dcount, zomega0);
	      err += nc_put_vara_double (dataFile, yomega,  dstart, dcount, zomega);
	      err += nc_put_vara_double (dataFile, yW,      dstart, dcount, zW);
	      err += nc_put_vara_double (dataFile, yWv,     dstart, dcount, zWv);
	      err += nc_put_vara_double (dataFile, yPsim,   dstart, dcount, zPsim);
	      err += nc_put_vara_double (dataFile, yPsip,   dstart, dcount, zPsip);
	      err += nc_put_vara_double (dataFile, yPsivm,  dstart, dcount, zPsivm);
	      err += nc_put_vara_double (dataFile, yPsivp,  dstart, dcount, zPsivp);
	      err += nc_put_vara_double (dataFile, ydP ,    dstart, dcount, zdP);
	      err += nc_put_vara_double (dataFile, ydP0,    dstart, dcount, zdP0);
	      err += nc_put_vara_double (dataFile, yTRmp,   dstart, dcount, zTRmp);
	      err += nc_put_vara_double (dataFile, yTWall,  dstart, dcount, zTWall);
	      err += nc_put_vara_double (dataFile, yTTear,  dstart, dcount, zTTear);

	      if (err != 0)
		{
		  printf ("PHASE::IslandDynamics: Error writing data to Outputs/Stage5.nc\n");
		  exit (1);
		}
	    }

	  if (cnt%10 == 0)
	    printf ("t(ms) = %10.3e h(ms) = %10.3e h/tau_A = %10.3e irmp(kA) = %10.3e prmp/pi = %10.3e\n", t*tau_A*1.e3, h*tau_A*1.e3, h, irmp, prmp /M_PI);
	  cnt++;

	  fflush (stdout);
	}
    }
  while (t < Tend);
  RK4Fixed (t, y, Tend - t);

  printf ("t(ms) = %10.3e h(ms) = %10.3e h/tau_A = %10.3e irmp(kA) = %10.3e prmp/pi = %10.3e\n", t*tau_A*1.e3, h*tau_A*1.e3, h, irmp, prmp /M_PI);

  err += nc_put_var_double (dataFile, ytlock,  ztlock);
  err += nc_put_var_double (dataFile, yIlock,  zIlock);
  err += nc_put_var_double (dataFile, ynelock, znelock);
  err += nc_put_var_double (dataFile, yTelock, zTelock);
  err += nc_put_var_double (dataFile, yCplock, zCplock);
  err += nc_put_var_double (dataFile, yB0lock, zB0lock);
  err += nc_put_var_double (dataFile, yR0lock, zR0lock);
  err += nc_put_var_double (dataFile, ywllock, zwllock);
  err += nc_put_var_double (dataFile, ywelock, zwelock);
  err += nc_put_var_double (dataFile, ywnlock, zwnlock);
  err += nc_put_var_double (dataFile, ybrlock, zbrlock);
  err += nc_put_var_double (dataFile, yWlock,  zWlock);

  if (err != 0)
    {
      printf ("PHASE::IslandDynamics: Error writing locking data to Outputs/Stage5.nc\n");
      exit (1);
    }
  
  err += nc_close (dataFile);

  if (err != 0)
    {
      printf ("PHASE::IslandDynamics: Error closing Outputs/Stage5.nc\n");
      exit (1);
    }

  // Save calculation for restart
  Save ();
 
  // Clean up
  delete[] zTime;   delete[] zIrmp;   delete[] zPrmp;   delete[] zphi;    delete[] zvph;    
  delete[] zomega0; delete[] zomega;  delete[] zPsim;   delete[] zPsip;   delete[] zW;      
  delete[] zPsivm;  delete[] zPsivp;  delete[] zWv;     delete[] zdP;     delete[] zdP0;    
  delete[] tstart;  delete[] tcount;  delete[] dstart;  delete[] dcount;  delete[] mk_x;
  delete[] PsiN_x;  delete[] zTRmp;   delete[] zTWall;  delete[] zTTear;  delete[] zIU;
  delete[] zIM;     delete[] zIL;     delete[] zPU;     delete[] zPM;     delete[] zPL;
  delete[] ztlock;  delete[] zIlock;  delete[] znelock; delete[] zTelock; delete[] zCplock;
  delete[] zB0lock; delete[] zR0lock; delete[] zwllock; delete[] zwelock; delete[] zwnlock;
  delete[] zbrlock; delete[] zWlock;
}

// ###################################
// Function to calculate irmp and prmp
// ###################################
void Phase::CalcRMP (double t)
{
  // Type 1 programmed waveform
  if (TYPE == 1)
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
		  irmp = ((t - TT (i)) * ICTRL[i+1] + (TT (i+1) - t) * ICTRL[i]) /(TT (i+1) - TT (i));
		  prmp = ((t - TT (i)) * PCTRL[i+1] + (TT (i+1) - t) * PCTRL[i]) /(TT (i+1) - TT (i));
		}
	    }
	}
    }

  // Type 2 spike waveform
  if (TYPE == 2)
    {
      if (t < SSTART)
	{
	  irmp = BACK;
	  prmp = RPHA; 
	}
      else if (t >= SSTART && t <= SEND)
	{
	  irmp = SAMP * cos ((t - SSTART) * WAMOD);
	  prmp = RPHA + (t - SSTART) * WPMOD;
	}
      else
	{
	  irmp = BACK;
	  prmp = RPHA;
	}
    }

  // Type 3 repeated rmap waveform
  if (TYPE == 3)
    {
      int    i    = int (t /RPERIOD);
      double toff = t = double (i) * RPERIOD;

      irmp = RSTART + (REND - RSTART) * toff /RPERIOD;
      prmp = RPHA;
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
      y (cnt) = Uk (j); cnt++;
      y (cnt) = Vk (j); cnt++;
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
      Xk    (j) = y (cnt); cnt++;
      Yk    (j) = y (cnt); cnt++;
      Uk    (j) = y (cnt); cnt++;
      Vk    (j) = y (cnt); cnt++;
      Psik  (j) = sqrt (Xk (j) * Xk (j) + Yk (j) * Yk (j));
      phik  (j) = atan2 (Yk (j), Xk (j));
      Psiwk (j) = sqrt (Uk (j) * Uk (j) + Vk (j) * Vk (j));
      phiwk (j) = atan2 (Vk (j), Uk (j));
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
		     Array<double,1> UkRHS, Array<double,1> VkRHS,
		     Array<double,2> alphakpRHS, Array<double,2> betakpRHS,
		     Array<double,1> dydt)
{
  int cnt = 0;
  for (int j = 0; j < nres; j++)
    {
      dydt (cnt) = XkRHS (j); cnt++;
      dydt (cnt) = YkRHS (j); cnt++;
      dydt (cnt) = UkRHS (j); cnt++;
      dydt (cnt) = VkRHS (j); cnt++;
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
	  double w = (0.8227/2.) * 4. * R_0 * fack (j) * sqrt (fabs (Psik (j)));

	  double facTe = (etae (j) /(1. + etae (j))) * (WcrTek (j)*WcrTek (j) /(WcrTek (j)*WcrTek (j) + w*w));
	  double facne = (1.       /(1. + etae (j))) * (Wcrnek (j)*Wcrnek (j) /(Wcrnek (j)*Wcrnek (j) + w*w));
	  double facTi = (etai (j) /(1. + etai (j))) * (w*w                   /(WcrTik (j)*WcrTik (j) + w*w));
	  double facni = (1.       /(1. + etai (j))) * (w*w                   /(Wcrnek (j)*Wcrnek (j) + w*w));
	  
	  om = wke (j) + (wkl (j) - wke (j)) * (facTe + facne) + (wkn (j) - wke (j)) * (facTi + facni);
	}
      else if (FREQ == 2)
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

// ##################################################################
// Function to calculate electromagnetic torques at rational surfaces
// ##################################################################
void Phase::GetElectromagneticTorques (Array<double,1 >& y, double* T_Rmp, double* T_Wall, double* T_Tear)
{
  Unpack (y);

  double normalization = - 2.*M_PI*M_PI * R_0*R_0*R_0 * B_0*B_0 * double(ntor(0)) /mu_0;
  for (int j = 0; j < nres; j++)
    {
      T_Rmp [j] = normalization * EEh (j, j) * Psik (j) * chi (j) * sin (phik (j) - zeta (j));
      T_Wall[j] = normalization * Sigmaw (j) * Psik (j) * Psiwk (j) * sin (phik (j) - phiwk (j));
      T_Tear[j] = 0.;
      for (int k = 0; k < nres; k++)
	T_Tear[j] +=  normalization * EEh (j, k) * Psik (j) * Psik (k) * sin (phik (j) - phik (k) - xih (j, k));
    }
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
      double W      = R_0 * GetIslandWidth (j) /dPsiNdr (j) /a (j)/R_0 /rk (j);
      double Wk     = 0.8227 * R_0 * GetIslandWidth (j) /dPsiNdr (j) /2. /a (j)/R_0 /rk (j);
      double WTek   = 0.8227 * WcrTek (j) /2. /a (j)/R_0 /rk (j);
      double WTik   = 0.8227 * WcrTik (j) /2. /a (j)/R_0 /rk (j);
      double Wnek   = 0.8227 * Wcrnek (j) /2. /a (j)/R_0 /rk (j);
      double rhek   = 0.8227 * rhothe (j) /2. /a (j)/R_0 /rk (j);
      double rhik   = 0.8227 * rhothi (j) /2. /a (j)/R_0 /rk (j);
      double facbTe = (etae (j) /(1. + etae (j))) * Wk /(WTek*WTek + rhek*rhek + Wk*Wk);
      double facbne = (1.       /(1. + etae (j))) * Wk /(Wnek*Wnek + rhek*rhek + Wk*Wk);
      double facbTi = (etai (j) /(1. + etai (j))) * Wk /(WTik*WTik + rhik*rhik + Wk*Wk);
      double facbni = (1.       /(1. + etai (j))) * Wk /(Wnek*Wnek + rhik*rhik + Wk*Wk);
      double faccTe = (etae (j) /(1. + etae (j))) * Wk /(WTek*WTek + Wk*Wk);
      double faccne = (1.       /(1. + etae (j))) * Wk /(Wnek*Wnek + Wk*Wk);
      double faccTi = (etai (j) /(1. + etai (j))) * Wk /(WTik*WTik + Wk*Wk);
      double faccni = (1.       /(1. + etai (j))) * Wk /(Wnek*Wnek + Wk*Wk);
      double facc   =   (nek (j) /(nek (j) + nik (j))) * (faccTe + faccne)
	              + (nik (j) /(nek (j) + nik (j))) * (faccTi + faccni);
      double facpTi = (etai (j) /(1. + etai (j))) * Wk /(WTik*WTik + Wk*Wk) /(WTik*WTik + Wk*Wk);
      double facpni = (1.       /(1. + etai (j))) * Wk /(Wnek*Wnek + Wk*Wk) /(Wnek*Wnek + Wk*Wk);
 
      double boot = alphabe (j) * (facbTe + facbne) + alphabi (j) * (facbTi + facbni);
      double curv = alphac  (j) * facc;
      double polz = alphap  (j) * (facpTi + facpni);
      double wall = Sigmaw (j) * Sigmaw (j) /Deltaw (j);
      double poem = ((Poem1 (j) * W * log (1./W) + Poem2 (j) * W) * WTek*WTek + Poem3 (j) * W * Wk*Wk) /(WTek*WTek + Wk*Wk);
      
      Cosk (j) = EEh (j, j) * chi (j) * cos (zeta (j)) + (boot + curv + polz + wall - poem) * Xk (j) + Sigmaw (j) * Uk (j);
      Sink (j) = EEh (j, j) * chi (j) * sin (zeta (j)) + (boot + curv + polz + wall - poem) * Yk (j) + Sigmaw (j) * Vk (j);
      sink (j) = EEh (j, j) * chi (j) * sin (phik (j) - zeta (j)) + Sigmaw (j) * Psiwk (j) * sin (phik (j) - phiwk (j));

      for (int k = 0; k < nres; k++)
	{
	  Cosk (j) += EEh (j, k) * (cos (xih (j, k)) * Xk (k) - sin (xih (j, k)) * Yk (k));
	  Sink (j) += EEh (j, k) * (cos (xih (j, k)) * Yk (k) + sin (xih (j, k)) * Xk (k));
	  sink (j) += EEh (j, k) * Psik (k) * sin (phik (j) - phik (k) - xih (j, k));
	}
    }

  Array<double,1> XkRHS      (nres);
  Array<double,1> YkRHS      (nres);
  Array<double,1> UkRHS      (nres);
  Array<double,1> VkRHS      (nres);
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

      UkRHS (j) = (Deltaw (j) * Uk (j) + Sigmaw (j) * Xk (j)) /tauw;
      VkRHS (j) = (Deltaw (j) * Vk (j) + Sigmaw (j) * Yk (j)) /tauw;
      
      for (int i = 0; i < NFLOW; i++)
	{
	  alphakpRHS (j, i) = (torp (j, i) * Psik (j) * sink (j)
			       - (j1p (i)*j1p (i) /taumk (j) + 1. /tautk (j) + 1. /tauxk (j)) * alphakp (j, i))
	    /(1. + 2.*qhatk (j)*qhatk (j));

	  betakpRHS (j, i) = tort (j, i) * Psik (j) * sink (j)
	    - (j0p (i)*j0p (i) /taumk (j)  + 1. /tauxk (j)) * betakp (j, i);
	}

      vp (j) = (YkRHS (j) * Xk (j) - XkRHS (j) * Yk (j)) / (Xk (j)*Xk (j) + Yk(j) * Yk(j));
    }

  PackRhs (XkRHS, YkRHS, UkRHS, VkRHS, alphakpRHS, betakpRHS, dydt);
}

// ######################################################################
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
// ######################################################################
void Phase::RK4Adaptive (double& x, Array<double,1>& y, double& h, double& t_err, double acc, double S, int& rept,
			 int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
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
      //printf ("Phase::RK4Adpative: Warning - |h| < hmin at x = %10.3e\n", x);
      //exit (1);
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %10.3e hin = %10.3e err = %10.3e acc = %10.3e hest = %10.3e hout = %10.3ey count = %3d\n", 
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
//  x ... independent variable
//  y ... array of dependent variables
//  h ... step-length
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
