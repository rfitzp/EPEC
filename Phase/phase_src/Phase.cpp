// Phase.cpp

// PROGRAM ORGANIZATION:
// 
//       Phase:: Phase          ()
// void  Phase:: Solve          (int _STAGE2, int _INTF, int _INTN, int _INTU, int _OLD, double _TIME)
// void  Phase:: Read_Data      (int _STAGE2, int _INTF, int _INTN, int _INTU, int _OLD, double _TIME)
// void  Phase:: Scan_Shift     ()
// void  Phase:: Calc_Velocity  ()
// void  Phase:: Initialize     ()
// void  Phase:: Save           ()
// void  Phase:: IslandDynamics ()
// void  Phase:: CalcRMP        (double t)
// void  Phase:: CalcChiZeta    (double t)
// void  Phase:: Pack           (Array<double,1> y)
// void  Phase:: Unpack         (Array<double,1> y)
// void  Phase:: PackRhs        (Array<double,1> PsikRHS, Array<double,1> phikRHS,
//		                 Array<double,2> alphakpRHS, Array<double,2> betakpRHS, Array<double,1> dydt)
// void  Phase:: Rhs            (double t, Array<double,1>& y, Array<double,1>& dydt)
// void  Phase:: RK4Adaptive    (double& x, Array<double,1>& y, double& h, 
//			         double& t_err, double acc, double S, int& rept,
//			         int maxrept, double h_min, double h_max, int flag, int diag, FILE* file)
// void  Phase:: RK4Fixed       (double& x, Array<double,1>& y, double h)
// FILE* Phase:: OpenFilew      (char* filename)
// FILE* Phase:: OpenFilew      (char* filename)
// FILE* Phase:: OpenFilea      (char* filename)

#include "Phase.h"

// ###########
// Constructor
// ###########
Phase::Phase ()
{
  // --------------------------------------------------
  // Set default values adaptive integration parameters
  // --------------------------------------------------
  h0       = 1.e-4;
  acc      = 1.e-11;
  hmin     = 1.e-10;
  hmax     = 1.e2;
  maxrept  = 50;

  // ----------------------
  // Set physical constants
  // ----------------------
  e         = GSL_CONST_MKSA_ELECTRON_CHARGE;
  epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
  mu_0      = GSL_CONST_MKSA_VACUUM_PERMEABILITY;
  m_p       = GSL_CONST_MKSA_MASS_PROTON;
  m_e       = GSL_CONST_MKSA_MASS_ELECTRON;
}

// #########################
// Function to solve problem
// #########################
void Phase::Solve (int _STAGE2, int _INTF, int _INTN, int _INTU, int _OLD, double _TIME)
{
  // Read data
  Read_Data (_STAGE2, _INTF, _INTN, _INTU, _OLD, _TIME);

  // Scan RMP phase shift
  Scan_Shift ();

  if (STAGE2 != 0)
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
  gsl_vector_complex_free (ChiU);    gsl_vector_complex_free (ChiL);
}

// #####################
// Function to read data
// #####################
void Phase::Read_Data (int _STAGE2, int _INTF, int _INTN, int _INTU, int _OLD, double _TIME)
{
  // ......................................
  // Set default values of input parameters
  // ......................................
  NFLOW  = 200;
  STAGE2 = 1;
  INTF   = 0;
  INTN   = 0;
  INTU   = 0;
  OLD    = 0;
  DT     = 1.e-5;
  TIME   = 0.;

  // Read data from namelist.txt
  printf ("...............................\n");
  printf ("Reading data from namelist.txt:\n");
  printf ("...............................\n");

  TCTRL = new double[MAXCONTROLPOINTNUMBER];
  ICTRL = new double[MAXCONTROLPOINTNUMBER];
  PCTRL = new double[MAXCONTROLPOINTNUMBER];

  NameListRead (&NFLOW, &STAGE2, &INTF, &INTN, &INTU, &OLD, &DT, &TIME, &NCTRL, TCTRL, ICTRL, PCTRL);

  // Override namelist values with command line options
  if (_STAGE2 > -1)
    STAGE2 = _STAGE2;
  if (_INTF > 0.)
    INTF = _INTF;
  if (_INTN > 0.)
    INTN = _INTN;
  if (_INTU > 0.)
    INTU = _INTU;
  if (_OLD > -1)
    OLD = _OLD;
  if (_TIME > 0.)
    TIME = _TIME;
 
  // .............................
  // Output calculation parameters
  // .............................
  printf ("NFLOW = %4d STAGE2 = %2d INTF = %2d INTN = %2d INTU = %2d OLD = %2d DT = %11.4e TIME = %11.4e NCTRL = %4d\n",
	  NFLOW, STAGE2, INTF, INTN, INTU, OLD, DT, TIME, NCTRL);

  FILE* monitor = OpenFilea ((char*) "../monitor.txt");
  fprintf (monitor, "NFLOW = %4d STAGE2 = %2d INTF = %2d INTN = %2d INTU = %2d OLD = %2d DT = %11.4e TIME = %11.4e NCTRL = %4d\n",
	   NFLOW, STAGE2, INTF, INTN, INTU, OLD, DT, TIME, NCTRL);
  fclose (monitor);
  
  // ............
  // Sanity check
  // ............
  if (NFLOW <= 0)
    {
      printf ("PHASE:: NFLOW must be positive\n");
      exit (1);
    }
  if (DT <= 0.)
    {
      printf ("PHASE:: DT must be positive\n");
      exit (1);
    }

  for (int i = 0; i < NCTRL; i++)
    printf ("T = %11.4e  IRMP = %11.4e  PRMP/pi = %11.4e\n", TCTRL[i], ICTRL[i], PCTRL[i]/M_PI);

  TT.resize (NCTRL);

  // .................
  // Interpolate fFile
  // .................
  if (INTF > 1 && TIME > 0.)
    {
      system ("rm -rf ../Flux/fFile");
      
      // Read fFile data
      char           filename[MAXFILENAMELENGTH];
      vector<string> fFileName;
      double         filetime;
      vector<double> fFileTime;
      int            fFileNumber = 0;
      
      printf ("Reading fFile data:\n");

      FILE* file = OpenFiler ((char*) "../Flux/fFileIndex");

      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  fFileName.push_back (filename);
	  fFileTime.push_back (filetime);
	}
      fFileNumber = fFileName.size ();

      fclose (file);
      
      // Interpolate fFiles
      fFileInterp (fFileName, fFileTime, fFileNumber, TIME);
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
 
  FILE* file = OpenFiler ((char*) "../Flux/fFile");
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
	      &R_0, &B_0, &inr, &q95, &r95, &qlim, &rlim, &q0, &qa, &NPSI, &ini, &nres) != 12)
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

  for (int j = 0; j < nres; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &ini, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr, &inr) != 11)
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
 
  EI = gsl_vector_complex_alloc (nres);
  EO = gsl_vector_complex_alloc (nres);
  
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

  epsi.resize (nres);
  sigi.resize (nres);
  epso.resize (nres);
  sigo.resize (nres);
  
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
  if (INTN > 1 && TIME > 0.)
    {
      system ("rm -rf ../Neoclassical/nFile");
      
      // Read nFile data
      char           filename[MAXFILENAMELENGTH];
      vector<string> nFileName;
      vector<double> nFileTime;
      double         filetime;
      int            nFileNumber = 0;
      
      printf ("Reading nFile data:\n");

      FILE* file = OpenFiler ((char*) "../Neoclassical/nFileIndex");
 
      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  nFileName.push_back (filename);
	  nFileTime.push_back (filetime);
	}
      nFileNumber = nFileName.size ();

      fclose (file);
 
      // Interpolate nFiles
      nFileInterp (nFileName, nFileTime, nFileNumber, TIME);
    }
  
  // ...................................
  // Read data from program NEOCLASSICAL
  // ...................................
  printf (".......................................\n");
  printf ("Reading data from program NEOCLASSICAL:\n");
  printf (".......................................\n");
 
  file = OpenFiler ((char*) "../Neoclassical/nFile");
  int nresn;
  if (fscanf (file, "%d %lf", &nresn, &tau_A) != 2)
    {
      printf ("PHASE::Error reading nFile (1)\n");
      exit (1);
    }
  if (nresn != nres)
    {
      printf ("PHASE:: Error - nres mismatch\n");
      exit (1);
    }

  printf ("tau_A = %11.4e\n", tau_A);

  for (int i = 0; i < NCTRL; i++)
    TT (i) = TCTRL[i] /tau_A;
  dTT = DT /tau_A;

  mk.resize     (nres); ntor.resize   (nres); rk.resize     (nres);
  qk.resize     (nres); rhok.resize   (nres); a.resize      (nres);
  Sk.resize     (nres); wk.resize     (nres); taumk.resize  (nres);
  tautk.resize  (nres); fack.resize   (nres); delk.resize   (nres);
  wkl.resize    (nres); wke.resize    (nres); dnedrk.resize (nres);
  dTedrk.resize (nres); Wcrnek.resize (nres); WcrTek.resize (nres);
  akk.resize    (nres); gk.resize     (nres); dPsidr.resize (nres);

  for (int j = 0; j < nres; j++)
    if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&mk(j), &ntor(j), &rk(j), &qk(j), &rhok(j), &a(j), &Sk(j),
		&wk(j), &taumk(j), &tautk(j), &fack(j), &delk(j), &wkl(j),
		&wke(j), &dnedrk(j), &dTedrk(j), &Wcrnek(j), &WcrTek(j), &akk(j), &gk(j), &dPsidr(j)) != 21)
      {
	printf ("PHASE::Error reading nFile (2)\n");
	exit (1);
      }
  fclose (file);

  for (int j = 0; j < nres; j++)
    printf ("m = %3d h_r = %10.3e q = %10.3e g = %10.3e akk = %10.3e h_rho = %10.3e h_a = %10.3e S = %10.3e h_w0 = %10.3e h_tauM = %10.3e h_tauth = %10.3e h_del = %10.3e\n",
	    mk(j), rk(j), qk(j), gk(j), akk(j), rhok(j), a(j), Sk(j), wk(j), taumk(j), tautk(j), delk(j));

  // .............................
  // Interpolate uFiles and lFiles
  // .............................
  if (INTU > 0  && TIME > 0.)
    {
      system ("rm -rf uFile lFile");

      // Read uFile data
      char           ufilename[MAXFILENAMELENGTH];
      vector<string> uFileName;
      double         ufiletime;
      vector<double> uFileTime;
      int            uFileNumber = 0;
      
      printf ("Reading uFile data:\n");

      file = OpenFiler ((char*) "uFileIndex");
 
      while (fscanf (file, "%s %lf %d %d %d", &ufilename, &ufiletime) == 2)
	{
	  uFileName.push_back (ufilename);
	  uFileTime.push_back (ufiletime);
	}
      uFileNumber = uFileName.size ();
      
      fclose (file);
 
      // Interpolate uFiles
      uFileInterp (uFileName, uFileTime, uFileNumber, TIME);

      // Read lFile data
      char           lfilename[MAXFILENAMELENGTH];
      vector<string> lFileName;
      double         lfiletime;
      vector<double> lFileTime;
      int            lFileNumber = 0;
      
      printf ("Reading lFile data:\n");

      file = OpenFiler ((char*) "lFileIndex");
      
      while (fscanf (file, "%s %lf", &lfilename, &lfiletime) == 2)
	{
	  lFileName.push_back (lfilename);
	  lFileTime.push_back (lfiletime);
	}
      lFileNumber = lFileName.size ();
      
      fclose (file);

      // Interpolate lFiles
      lFileInterp (lFileName, lFileTime, lFileNumber, TIME);
    }

  // ...........................
  // Read data from program GPEC
  // ...........................
  printf ("...............................\n");
  printf ("Reading data from program GPEC:\n");
  printf ("...............................\n");
  
  char line[MAXULFILELINELENGTH]; char line1[MAXULFILELINELENGTH];
  double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;
  double* QIN = new double[nres];
  double* DRE = new double[nres];
  double* DIM = new double[nres];
  double* CRE = new double[nres];
  double* CIM = new double[nres];
  double* WWW = new double[nres];
  DeltaU = gsl_vector_complex_calloc (nres);
  DeltaL = gsl_vector_complex_calloc (nres);
  ChiU   = gsl_vector_complex_calloc (nres);
  ChiL   = gsl_vector_complex_calloc (nres);

  double SCALEFACTOR = 2.;
  printf ("SCALEFACTOR = %11.4e\n", SCALEFACTOR);
  
  printf ("Upper coil:\n");
  file = OpenFiler ((char*) "uFile");

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
      printf ("PHASE:: Warning - nsing < nres\n");
      nres1 = nsingu;
    }
  if (nsingu > nres)
    {
      printf ("PHASE:: Warning - nsing > nres\n");
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
	DRE[i] = v9;
	DIM[i] = v10;
	WWW[i] = 2.*v11;

	DRE[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;
	DIM[i] /= 2.*M_PI * qk(i) /SCALEFACTOR/SCALEFACTOR;

	CRE[i] =   DIM[i] * (rk(i) * a(i)) * (rk(i) * a(i)) * gk(i) /double (mk(i)) /(akk(i) + rk(i) * rk(i) * a(i) * a(i) /qk(i) /qk(i)) /EEh(i, i);
	CIM[i] = - DRE[i] * (rk(i) * a(i)) * (rk(i) * a(i)) * gk(i) /double (mk(i)) /(akk(i) + rk(i) * rk(i) * a(i) * a(i) /qk(i) /qk(i)) /EEh(i, i);

	gsl_vector_complex_set (DeltaU, i, gsl_complex_rect (DRE[i], DIM[i]));
	gsl_vector_complex_set (ChiU,   i, gsl_complex_rect (CRE[i], CIM[i]));

	double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiU, i));
	double WUNRE = 4. * fack(i) * sqrt (Psi) * dPsidr (i);
	double WFULL = sqrt (FFh(i, i) * EEh(i, i)) * WUNRE;

	printf ("q = %11.4e  Delta = (%11.4e, %11.4e)  Chi = (%11.4e, %11.4e)  W_UNRE = %11.4e  W_UNRE/W_GPEC = %11.4e  W_FULL/W_GPEC = %11.4e\n",
		QIN[i], DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);

      }
  fclose (file);

  printf ("Lower coil:\n");
  file = OpenFiler ((char*) "lFile");
  
  for (int i = 0; i < 5; i++)
    fgets (line, MAXULFILELINELENGTH, file);
  fgets (line1, MAXULFILELINELENGTH, file);
  for (int i = 0; i < 2; i++)
    fgets (line, MAXULFILELINELENGTH, file);

  token = strtok (line1, " "); 
  token = strtok (NULL, " "); 
  token = strtok (NULL, " ");
  int nsingl = atoi (token);

  if (nsingl != nsingu)
    {
      printf ("PHASE:: Error - nsingl != nsingu\n");
      exit (1);
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
	DRE[i] = v9;
	DIM[i] = v10;
	WWW[i] = 2.*v11;

	DRE[i] /= 2.*M_PI * qk(i) / SCALEFACTOR/SCALEFACTOR;
	DIM[i] /= 2.*M_PI * qk(i) / SCALEFACTOR/SCALEFACTOR;

	CRE[i] =   DIM[i] * (rk(i) * a(i)) * (rk(i) * a(i)) * gk(i) /double (mk(i)) /(akk(i) + rk(i) * rk(i) * a(i) * a(i) /qk(i) /qk(i)) /EEh(i, i);
	CIM[i] = - DRE[i] * (rk(i) * a(i)) * (rk(i) * a(i)) * gk(i) /double (mk(i)) /(akk(i) + rk(i) * rk(i) * a(i) * a(i) /qk(i) /qk(i)) /EEh(i, i);

	gsl_vector_complex_set (DeltaL, i, gsl_complex_rect (DRE[i], DIM[i]));
	gsl_vector_complex_set (ChiL,   i, gsl_complex_rect (CRE[i], CIM[i]));

	double Psi   = gsl_complex_abs (gsl_vector_complex_get (ChiL, i));
	double WUNRE = 4. * fack(i) * sqrt (Psi) * dPsidr(i);
	double WFULL = sqrt (FFh(i, i) * EEh(i, i)) * WUNRE;

	printf ("q = %11.4e  Delta = (%11.4e, %11.4e)  Chi = (%11.4e, %11.4e)  W_UNRE = %11.4e  W_UNRE/W_GPEC = %11.4e  W_FULL/W_GPEC = %11.4e\n",
		QIN[i], DRE[i], DIM[i], CRE[i], CIM[i], WUNRE, WUNRE/WWW[i], WFULL/WWW[i]);
      }
  fclose (file);

  file = OpenFilea ((char*) "../IslandDynamics/Stage1/q.txt");
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", q0, q95, qa, qlim, TIME);
  fclose (file);

  file = OpenFilea ((char*) "../IslandDynamics/Stage1/omega0.txt");
  for (int j = 0; j < nres; j++)
     fprintf (file, "%3d %16.9e %16.9e %16.9e %16.9e %16.9e\n", mk(j), rk(j), wkl(j) /tau_A/1.e3, wk(j) /tau_A/1.e3, wke(j) /tau_A/1.e3, TIME);
  fclose (file);

  delete[] QIN; delete[] DRE; delete[] DIM; delete[] CRE; delete[] CIM; delete[] WWW;
 }

// #####################################################################
// Function to calculate vacuum flux versus upper/lower coil phase-shift
// #####################################################################
void Phase::Scan_Shift ()
{
  FILE* file1 = OpenFilew ((char*) "Stage1/chi.txt");
  FILE* file2 = OpenFilew ((char*) "Stage1/zeta.txt");
  FILE* file3 = OpenFilew ((char*) "Stage1/vac.txt");

  int    I   = 360;
  double one = 1.;
  for (int i = 0; i <= I; i++)
    {
      double      pha = double (i) * 4.*M_PI /double (I);
      gsl_complex eik = gsl_complex_polar (one, pha);

      fprintf (file1, "%e", pha/M_PI);
      fprintf (file2, "%e", pha/M_PI);
      fprintf (file3, "%e", pha/M_PI);
      for (int j = 0; j < nres; j++)
	{
	  gsl_complex hu = gsl_vector_complex_get (ChiU, j);
	  gsl_complex hl = gsl_vector_complex_get (ChiL, j);
	  hl             = gsl_complex_mul (hl, eik);
	  gsl_complex h  = gsl_complex_add (hu, hl);

	  double chi     =   gsl_complex_abs (h);
	  double zeta    = - gsl_complex_arg (h);
	  double wv      = 4. * R_0 * fack (j) * sqrt (chi);
	    
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

  FILE* file4 = OpenFilea ((char*) "../IslandDynamics/Stage1/vac.txt");

  for (int j = 0; j < nres; j++)
    {
      gsl_complex hu = gsl_vector_complex_get (ChiU, j);
      gsl_complex hl = gsl_vector_complex_get (ChiL, j);
      gsl_complex h  = gsl_complex_add (hu, hl);
      
      double chi     =   gsl_complex_abs (h);
      double zeta    = - gsl_complex_arg (h);
      double wv      = 4. * fack (j) * sqrt (chi) /a (j);
      
      fprintf (file4, "%3d %16.9e %16.9e %16.9e\n", mk (j), rk (j), wv, TIME);
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

      for (int i = 0; i < NFLOW; i++)
	{
	  alphakp (j, i) = 0.;
	  betakp  (j, i) = 0.;
	}
    }

  for (int j = 0; j < nres; j++)
    {
      lock (j) = 0;

      double sum = wk (j);
	  
      for (int k = 0; k < nres; k++)
	for (int i = 0; i < NFLOW; i++)
	  sum -= alphakp (k, i) * natp (j, k, i) + betakp (k, i) * natt (j, k, i);
      
      ww (j) = sum;
    }

  // Reinitialize from previous run
  if (OLD == 1)
    {
      FILE* file = OpenFiler ((char*) "sFile");

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
	if (fscanf (file, "%d %lf %lf %d %lf\n", &in, &_Psik(j), &_phik(j), &_lock(j), &_ww(j)) != 5)
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
	  Psik(j) = _Psik(j);
	  phik(j) = _phik(j);
	  lock(j) = _lock(j);
	  ww  (j) = _ww  (j);
	}

      for (int j = 0; j < nres; j++)
	for (int i = 0; i < NFLOW; i++)
	  {
	    alphakp (j, i) = _alphakp (j, i);
	    betakp  (j, i) = _betakp  (j, i);
	  }
    }
}

// ############################################
// Function to save island dynamics calculation
// ############################################
void Phase::Save ()
{
  FILE* file = OpenFilew ("sFile");

  fprintf (file, "%d %d\n", nres, NFLOW);

  for (int j = 0; j < nres; j++)
    fprintf (file, "%3d %19.6e %19.6e %2d %19.6e\n", j, Psik(j), phik(j), lock(j), ww(j));

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
  
  t = 0.;
  h = h0;
  Pack (y);

  int    cnt    = 0;
  double dt     = 0.; 
  FILE*  file   = OpenFilew ((char*) "Stage2/ntor.txt");
  FILE*  file1  = OpenFilew ((char*) "Stage2/Psi.txt");
  FILE*  file2  = OpenFilew ((char*) "Stage2/phi.txt");
  FILE*  file3  = OpenFilew ((char*) "Stage2/W.txt");
  FILE*  file4  = OpenFilew ((char*) "Stage2/omega.txt");
  FILE*  file5  = OpenFilew ((char*) "Stage2/rmp.txt");
  FILE*  file6  = OpenFilew ((char*) "Stage2/mirnov.txt");
  FILE*  file8  = OpenFilea ((char*) "../IslandDynamics/Stage1/scan.txt");  
  FILE*  file9  = OpenFilew ((char*) "Stage2/deltane.txt");
  FILE*  file10 = OpenFilew ((char*) "Stage2/deltaTe.txt");
  FILE*  file11 = OpenFilew ((char*) "Stage2/dne.txt");
  FILE*  file12 = OpenFilew ((char*) "Stage2/dTe.txt");
  FILE*  file13 = OpenFilew ((char*) "Stage2/rkminus.txt");
  FILE*  file14 = OpenFilew ((char*) "Stage2/rkplus.txt");
  FILE*  file15 = OpenFilew ((char*) "Stage2/rkminusne.txt");
  FILE*  file16 = OpenFilew ((char*) "Stage2/rkplusne.txt");
  FILE*  file17 = OpenFilew ((char*) "Stage2/rkminusTe.txt");
  FILE*  file18 = OpenFilew ((char*) "Stage2/rkplusTe.txt");

  fprintf (file, "%d\n", ntor (0));
  fclose (file);
  
  printf ("......................\n");
  printf ("Performing simulation:\n");
  printf ("......................\n");
  do
    {
      for (int j = 0; j < nres; j++)
	wwo (j) = ww (j);
      
      RK4Adaptive (t, y, h, t_err, acc, 2., rept, maxrept, hmin, hmax, 2, 0, NULL);
      dt += h;

      Unpack (y);

      for (int j = 0; j < nres; j++)
	{
	  double sum = wk (j);
	  
	  for (int k = 0; k < nres; k++)
	    for (int i = 0; i < NFLOW; i++)
	      sum -= alphakp (k, i) * natp (j, k, i) + betakp (k, i) * natt (j, k, i);
	  
	  ww (j) = sum;
	}

      for (int j = 0; j < nres; j++)
	if (ww (j) * wwo (j) < 0. && lock (j) == 0)
	  {
	    CalcRMP (t); 

	    printf ("m = %3d locks at t = %11.4e s  irmp = %11.4e kA  prmp/pi = %11.4e\n",
		    mk (j), t*tau_A, irmp, prmp /M_PI);
	    fprintf (file8, "%16.9e %3d %16.9e %16.9e %16.9e %16.9e %3d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
		     q95, mk (j), rk (j), t*tau_A, irmp, prmp /M_PI, nres, q0, qa,
		     wk(j)/tau_A/1.e3, wkl(j)/tau_A/1.e3, wke(j)/tau_A/1.e3, TIME);

	    lock (j) = 1;
	  }
      	  
      if (dt > dTT)
	{
	  dt = 0.;

	  fprintf (file1, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file1, "%16.9e ", Psik (j));
	  fprintf (file1, "\n");

	  fprintf (file2, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file2, "%16.9e ", atan2 (sin (phik (j)), cos (phik (j))) /M_PI);
	  fprintf (file2, "\n");

	  fprintf (file3, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file3, "%16.9e ", 4. * R_0 * fack (j) * sqrt (fabs (Psik (j))));
	  fprintf (file3, "\n");

	  fprintf (file4, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    fprintf (file4, "%16.9e ", ww (j) /tau_A/1.e3);
	  fprintf (file4, "\n");

	  CalcRMP (t);
	  fprintf (file5, "%16.9e %16.9e %16.9e\n", t*tau_A, irmp, prmp /M_PI);

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
	  
	  fprintf (file6, "%16.9e %16.9e %16.9e %16.9ee %16.9e\n", t*tau_A, sumci*B_0*1.e4, sumsi*B_0*1.e4, sumco*B_0*1.e4, sumso*B_0*1.e4);

	  fprintf (file9,  "%16.9e ", t*tau_A); fprintf (file10, "%16.9e ", t*tau_A); fprintf (file11, "%16.9e ", t*tau_A);
	  fprintf (file12, "%16.9e ", t*tau_A); fprintf (file13, "%16.9e ", t*tau_A); fprintf (file14, "%16.9e ", t*tau_A);
	  fprintf (file15, "%16.9e ", t*tau_A); fprintf (file16, "%16.9e ", t*tau_A); fprintf (file17, "%16.9e ", t*tau_A);
	  fprintf (file18, "%16.9e ", t*tau_A);
	  for (int j = 0; j < nres; j++)
	    {
	      double Wk       = 4. * R_0 * fack (j) * sqrt (fabs (Psik (j)));
	      double deltanek = (2./M_PI) * Wk *Wk*Wk /(Wk*Wk + Wcrnek(j) * Wcrnek(j));
	      double deltaTek = (2./M_PI) * Wk *Wk*Wk /(Wk*Wk + WcrTek(j) * WcrTek(j));
	      double dnek     = dnedrk (j) * deltanek;
	      double dTek     = dTedrk (j) * deltaTek;
	      fprintf (file9,  "%16.9e ", deltanek);
	      fprintf (file10, "%16.9e ", deltaTek);
	      fprintf (file11, "%16.9e ", dnek);
	      fprintf (file12, "%16.9e ", dTek);
	      fprintf (file13, "%16.9e ", rk(j) - Wk /2./a(j)/R_0);
	      fprintf (file14, "%16.9e ", rk(j) + Wk /2./a(j)/R_0);
	      fprintf (file15, "%16.9e ", rk(j) - deltanek /2./a(j)/R_0);
	      fprintf (file16, "%16.9e ", rk(j) + deltanek /2./a(j)/R_0);
	      fprintf (file17, "%16.9e ", rk(j) - deltaTek /2./a(j)/R_0);
	      fprintf (file18, "%16.9e ", rk(j) + deltaTek /2./a(j)/R_0);
	    }
	  fprintf (file9,  "\n"); fprintf (file10, "\n"); fprintf (file11, "\n"); fprintf (file12, "\n");
	  fprintf (file13, "\n"); fprintf (file14, "\n"); fprintf (file15, "\n"); fprintf (file16, "\n");
	  fprintf (file17, "\n"); fprintf (file18, "\n");

	  if (cnt%100 == 0)
	    printf ("t = %11.4e  irmp = %11.4e  prmp/pi = %11.4e\n", t*tau_A, irmp, prmp /M_PI);
	  cnt++;
	  
	  fflush (file1);  fflush (file2);  fflush (file3);
	  fflush (file4);  fflush (file5);  fflush (file6);
	  fflush (file8);  fflush (file9);  fflush (file10);
	  fflush (file11); fflush (file12); fflush (file13);
	  fflush (file14); fflush (file15); fflush (file16);
	  fflush (file17); fflush (file18);
	}
    }
  while (t < TT (NCTRL-1));
  fclose (file1);  fclose (file2);  fclose (file3);
  fclose (file4);  fclose (file5);  fclose (file6);
  fclose (file8);  fclose (file9);  fclose (file10);
  fclose (file11); fclose (file12); fclose (file13);
  fclose (file14); fclose (file15); fclose (file17);
  fclose (file17); fclose (file18);

  // Save calculation
  Save ();
  
  FILE* filew = OpenFilea ((char*) "../IslandDynamics/Stage1//omega.txt");
  for (int j = 0; j < nres; j++)
    fprintf (file, "%3d %16.9e %16.9e %16.9e %16.9e\n", mk (j), rk (j), wk (j) /tau_A/1.e3, ww (j) /tau_A/1.e3, TIME);
  fclose (file);
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
	      irmp = ((t - TT (i)) * ICTRL[i+1] + (TT (i+1) - t) * ICTRL[i]) /(TT (i+1) - TT (i));
	      prmp = ((t - TT (i)) * PCTRL[i+1] + (TT (i+1) - t) * PCTRL[i]) /(TT (i+1) - TT (i));
	    }
	}
    }
}

// ##########################################
// Function to calculate chi and zeta vectors
// ##########################################
void Phase::CalcChiZeta (double t)
{
  CalcRMP (t);
  
  double      one = 1.;
  gsl_complex eik = gsl_complex_polar (one, prmp);
  
  for (int j = 0; j < nres; j++)
    {
      gsl_complex hu = gsl_vector_complex_get (ChiU, j);
      gsl_complex hl = gsl_vector_complex_get (ChiL, j);
      hl             = gsl_complex_mul (hl, eik);
      gsl_complex h  = gsl_complex_add (hu, hl);

      chi  (j) =   gsl_complex_abs (h) * irmp;
      zeta (j) = - gsl_complex_arg (h);
    }
}

// ########################################################
// Function to pack simulation variables into single vector
// ########################################################
void Phase::Pack (Array<double,1> y)
{
  int cnt = 0;
  for (int j = 0; j  < nres; j++)
    {
      y (cnt) = Psik (j); cnt++;
      y (cnt) = phik (j); cnt++;
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
  for (int j = 0; j  < nres; j++)
    {
      Psik (j) = y (cnt); cnt++;
      phik (j) = y (cnt); cnt++;
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
void Phase::PackRhs (Array<double,1> PsikRHS, Array<double,1> phikRHS,
		     Array<double,2> alphakpRHS, Array<double,2> betakpRHS,
		     Array<double,1> dydt)
{
  int cnt = 0;
  for (int j = 0; j  < nres; j++)
    {
      dydt (cnt) = PsikRHS (j); cnt++;
      dydt (cnt) = phikRHS (j); cnt++;
      for (int i = 0; i < NFLOW; i++)  
	{
	  dydt (cnt) = alphakpRHS (j, i); cnt++;
	  dydt (cnt) = betakpRHS  (j, i); cnt++;
	}
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

  for (int j = 0; j < nres; j++)
    {
      Cosk (j) = EEh (j, j) * chi (j) * cos (phik (j) - zeta (j));
      Sink (j) = EEh (j, j) * chi (j) * sin (phik (j) - zeta (j));

      for (int k = 0; k < nres; k++)
	{
	  Cosk (j) += EEh (j, k) * Psik (k) * cos (phik (j) - phik (k) - xih (j, k));
	  Sink (j) += EEh (j, k) * Psik (k) * sin (phik (j) - phik (k) - xih (j, k));
	}
    }

  Array<double,1> PsikRHS    (nres);
  Array<double,1> phikRHS    (nres);
  Array<double,2> alphakpRHS (nres, NFLOW);
  Array<double,2> betakpRHS  (nres, NFLOW);

  for (int j = 0; j < nres; j++)
    {
      PsikRHS (j) = Cosk (j) /Sk (j) /(sqrt (fabs (Psik (j))) + delk (j));

      double sum = wk (j);
      for (int k = 0; k < nres; k++)
	for (int i = 0; i < NFLOW; i++)
	  sum -= alphakp (k, i) * natp (j, k, i) + betakp (k, i) * natt (j, k, i);
      phikRHS (j) = sum;
      
      for (int i = 0; i < NFLOW; i++)
	{
	  alphakpRHS (j, i) = (torp (j, i) * Psik (j) * Sink (j)
			       - (j1p (i)*j1p (i) /taumk (j) + 1. /tautk (j)) * alphakp (j, i))
	    /(1. + 2.*qk (j)*qk (j));

	  betakpRHS (j, i) = tort (j, i) * Psik (j) * Sink (j)
	    - (j0p (i)*j0p (i) /taumk (j)) * betakp (j, i);
	}
    }

  PackRhs (PsikRHS, phikRHS, alphakpRHS, betakpRHS, dydt);
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
          err2  =  fabs (y (i) - y1 (i));
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
  else if (h_est / h < 1. /S)
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
