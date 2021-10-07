// Rescale.cpp

#include "Rescale.h"

// ###########
// Constructor
// ###########
Rescale::Rescale ()
{
}

// ###############################
// Function to rescale equilibrium
// ###############################
void Rescale::RescaleEquilibrium ()
{
  // Read namelist
  NameListRead (&TYPE, &SCALE, &PSHIFT, &WSHIFT, &Q95, &OPOINT, &XPOINT);

  if (TYPE == 1)
    {
      // Rescale gFile
      double an  = SCALE;
      double anh = sqrt (an);
      gFileRescaleTypeI (&anh, &OPOINT, &XPOINT);

      printf ("\nType I gFile rescaling:  scale factor = %11.4e\n\n", anh);

      // Rescale pFile
      pFileRescaleType1 (an);
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (1.);
      printf ("\n");
    }
  else if (TYPE == 2)
    {
      // Rescale gFile
      double at  = SCALE;
      double ath = sqrt (at);
      gFileRescaleTypeI (&ath, &OPOINT, &XPOINT);

      printf ("\nType I gFile rescaling:  scale factor = %11.4e\n\n", ath);

      // Rescale pFile
      pFileRescaleType2 (at);
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (1.);
      printf ("\n");
    }
  else if (TYPE == 3)
    {
      // Rescale gFile
      double ar = SCALE;
      gFileRescaleTypeII (&ar, &OPOINT, &XPOINT);

      printf ("\nType II gFile rescaling:  scale factor = %11.4e\n\n", ar);

      // Rescale pFile
      pFileRescaleType3 (ar);
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (1.);
      printf ("\n");
    }
  else if (TYPE == 4)
    {
      // Rescale gFile
      double ac  = SCALE;
      double ach = 1.;
      gFileRescaleType0 ();

      printf ("\nType O gFile rescaling:  scale factor = %11.4e\n\n", ach);

      // Rescale pFile
      pFileRescaleType4 (ach);
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (ac);
      printf ("\n");
    }
  else if (TYPE == 5)
    {
      // Rescale gFile
      gFileRescaleTypeIII (&PSHIFT);

      printf ("\nType III gFile rescaling:  pressure shift = %11.4e\n\n", PSHIFT);

      // Rescale pFile
      {
	Field ine, ini;
	GetnFields (ine, ini);
	pFileRescaleType5 (PSHIFT, ine, ini);
      }
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (1.);
      printf ("\n");
      fflush (stdout);
    }
  else if (TYPE == 6)
    {
      // Rescale gFile
      gFileRescaleType0 ();

      double a = 1.;
      printf ("\nType 0 gFile rescaling:  scale factor = %11.4e\n\n", a);

      // Rescale pFile
      {
	Field RBp, R;
	GetRFields (RBp, R);
	pFileRescaleType6 (WSHIFT, RBp, R);
      }
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (1.);
      printf ("\n");
    }
  else if (TYPE == 7)
    {
      // Rescale gFile
      double q95_old, a1;   
      gFileRescaleType7 (&Q95, &OPOINT, &XPOINT, &q95_old, &a1);
      
      printf ("\nType IV/I gFile rescaling: q95_old = %11.4e  q95_new = %11.4e  pressure rescale factor = %11.4e\n\n", q95_old, Q95, a1);
      
      // Rescale pFile
      pFileRescaleType7 (a1);
      printf ("\n");

      // Rescale cFile
      cFileRescaleType4 (1.);
      printf ("\n");
    }
}

// #######################################################
// Function to extract 1/2/ne and 1/2/ni fields from pFile
// #######################################################
void Rescale::GetnFields (Field& ine, Field& ini)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field, ne, ni;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::GetFields: Error opening input pFile\n");
      exit (1);
    }

  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::GetnFields: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::GetnFields: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      if (strstr (s, "ne") != NULL)
	field.Copy (ne);
      else if (strstr (s, "ni") != NULL)
	field.Copy (ni);
    }
  while (1);

  fclose (file);

  double fac = 1.e20 * 1.60217662e-19;
  n = ne.GetN();
  ine.resize (n);
  ini.resize (n);
  for (int i = 0; i < n; i++)
    {
      double x1, y1, dydx1;
      ne.PullData (i, x1, y1, dydx1);
      double x2, y2, dydx2;
      ni.PullData (i, x2, y2, dydx2);
 
      double x3    = x1;
      double y3    = 1./2. /y1 /fac;
      double dydx3 = 0.;
      ine.PushData (i, x3, y3, dydx3);

      double x4    = x1;
      double y4    = 1./2. /y2/ fac;
      double dydx4 = 0.;
      ini.PushData (i, x4, y4, dydx4);
    }
}

// ###############################################
// Function to extract RBp and R fields from pFile
// ###############################################
void Rescale::GetRFields (Field& RBp, Field& R)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field, omeg, omgeb, er, vtor1;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::GetRFields: Error opening input pFile\n");
      exit (1);
    }

  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::GetRFields: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::GetRFields: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      if (strstr (s, "omeg") != NULL)
	field.Copy (omeg);
      else if (strstr (s, "omgeb") != NULL)
	field.Copy (omgeb);
      else if (strstr (s, "er") != NULL)
	field.Copy (er);
      else if (strstr (s, "vtor1") != NULL)
	field.Copy (vtor1);
    }
  while (1);

  fclose (file);

  n = omeg.GetN();
  RBp.resize (n);
  R.resize   (n);
  for (int i = 0; i < n; i++)
    {
      double x1, y1, dydx1;
      omeg.PullData (i, x1, y1, dydx1);
      double x2, y2, dydx2;
      omgeb.PullData (i, x2, y2, dydx2);
      double x3, y3, dydx3;
      er.PullData (i, x3, y3, dydx3);
      double x4, y4, dydx4;
      vtor1.PullData (i, x4, y4, dydx4);

      double x5    = x1;
      double y5    = y3 /y2;
      double dydx5 = 0.;
      RBp.PushData (i, x5, y5, dydx5);

      double x6    = x1;
      double y6    = y4 /y1;
      double dydx6 = 0.;
      R.PushData (i, x6, y6, dydx6);
    }
}

// #############################################
// Function to perform Type 1 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType1 (double an)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType1: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType1: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType1: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double An = 1.;
      if (strstr (s, "ne") != NULL)
	An = an;
      else if (strstr (s, "te") != NULL)
	An = 1.;
      else if (strstr (s, "ni") != NULL)
	An = an;
      else if (strstr (s, "ti") != NULL)
	An = 1.;
      else if (strstr (s, "nb") != NULL)
	An = an;
      else if (strstr (s, "pb") != NULL)
	An = an;
      else if (strstr (s, "ptot") != NULL)
	An = an;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "omegp") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "omgvb") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "omgpp") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "omgeb") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "er") != NULL)
	An = 1.;
      else if (strstr (s, "ommvb") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "ommpp") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "omevb") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "omepp") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "kpol") != NULL)
	An = 1./an;
      else if (strstr (s, "omghb") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "nz1") != NULL)
	An = an;
      else if (strstr (s, "vtor1") != NULL)
	An = 1./sqrt(an);
      else if (strstr (s, "vpol1") != NULL)
	An = 1./sqrt(an);

      field.Rescale (An);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));
      
      printf ("pFile field: %-34s:  rescale factor = %11.4e\n", s, An);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 2 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType2 (double at)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType2: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType2: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType2: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double At = 1.;
      if (strstr (s, "ne") != NULL)
	At = 1.;
      else if (strstr (s, "te") != NULL)
	At = at;
      else if (strstr (s, "ni") != NULL)
	At = 1.;
      else if (strstr (s, "ti") != NULL)
	At = at;
      else if (strstr (s, "nb") != NULL)
	At = 1.;
      else if (strstr (s, "pb") != NULL)
	At = at;
      else if (strstr (s, "ptot") != NULL)
	At = at;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	At = sqrt(at);
      else if (strstr (s, "omegp") != NULL)
	At = sqrt(at);
      else if (strstr (s, "omgvb") != NULL)
	At = sqrt(at);
      else if (strstr (s, "omgpp") != NULL)
	At = sqrt(at);
      else if (strstr (s, "omgeb") != NULL)
	At = sqrt(at);
      else if (strstr (s, "er") != NULL)
	At = at;
      else if (strstr (s, "ommvb") != NULL)
	At = sqrt(at);
      else if (strstr (s, "ommpp") != NULL)
	At = sqrt(at);
      else if (strstr (s, "omevb") != NULL)
	At = sqrt(at);
      else if (strstr (s, "omepp") != NULL)
	At = 1./sqrt(at);
      else if (strstr (s, "kpol") != NULL)
	At = 1.;
      else if (strstr (s, "omghb") != NULL)
	At = sqrt(at);
      else if (strstr (s, "nz1") != NULL)
	At = 1.;
      else if (strstr (s, "vtor1") != NULL)
	At = sqrt(at);
      else if (strstr (s, "vpol1") != NULL)
	At = sqrt(at);

      field.Rescale (At);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));
 
      printf ("pFile field: %-34s:  rescale factor = %11.4e\n", s, At);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 3 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType3 (double ar)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType3: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType3: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType3: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double Ar = 1.;
      if (strstr (s, "ne") != NULL)
	Ar = 1.;
      else if (strstr (s, "te") != NULL)
	Ar = 1.;
      else if (strstr (s, "ni") != NULL)
	Ar = 1.;
      else if (strstr (s, "ti") != NULL)
	Ar = 1.;
      else if (strstr (s, "nb") != NULL)
	Ar = 1.;
      else if (strstr (s, "pb") != NULL)
	Ar = 1.;
      else if (strstr (s, "ptot") != NULL)
	Ar = 1.;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "omegp") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "omgvb") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "omgpp") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "omgeb") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "er") != NULL)
	Ar = 1./ar;
      else if (strstr (s, "ommvb") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "ommpp") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "omevb") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "omepp") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "kpol") != NULL)
	Ar = 1./ar;
      else if (strstr (s, "omghb") != NULL)
	Ar = 1./ar/ar;
      else if (strstr (s, "nz1") != NULL)
	Ar = 1.;
      else if (strstr (s, "vtor1") != NULL)
	Ar = 1./ar;
      else if (strstr (s, "vpol1") != NULL)
	Ar = 1./ar;

      field.Rescale (Ar);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));
      
      printf ("pFile field: %-34s:  rescale factor = %11.4e\n", s, Ar);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 4 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType4 (double ac)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType4: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType4: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType4: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double Ac = 1.;
      if (strstr (s, "ne") != NULL)
	Ac = 1.;
      else if (strstr (s, "te") != NULL)
	Ac = 1.;
      else if (strstr (s, "ni") != NULL)
	Ac = 1.;
      else if (strstr (s, "ti") != NULL)
	Ac = 1.;
      else if (strstr (s, "nb") != NULL)
	Ac = 1.;
      else if (strstr (s, "pb") != NULL)
	Ac = 1.;
      else if (strstr (s, "ptot") != NULL)
	Ac = 1.;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	Ac = 1.;
      else if (strstr (s, "omegp") != NULL)
	Ac = 1.;
      else if (strstr (s, "omgvb") != NULL)
	Ac = 1.;
      else if (strstr (s, "omgpp") != NULL)
	Ac = 1.;
      else if (strstr (s, "omgeb") != NULL)
	Ac = 1.;
      else if (strstr (s, "er") != NULL)
	Ac = 1.;
      else if (strstr (s, "ommvb") != NULL)
	Ac = 1.;
      else if (strstr (s, "ommpp") != NULL)
	Ac = 1.;
      else if (strstr (s, "omevb") != NULL)
	Ac = 1.;
      else if (strstr (s, "omepp") != NULL)
	Ac = 1.;
      else if (strstr (s, "kpol") != NULL)
	Ac = 1.;
      else if (strstr (s, "omghb") != NULL)
	Ac = 1.;
      else if (strstr (s, "nz1") != NULL)
	Ac = 1.;
      else if (strstr (s, "vtor1") != NULL)
	Ac = 1.;
      else if (strstr (s, "vpol1") != NULL)
	Ac = 1.;

      field.Rescale (Ac);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));
 
      printf ("pFile field: %-34s:  rescale factor = %11.4e\n", s, Ac);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 5 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType5 (double ap, Field& ine, Field& ini)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType5: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType5: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType5: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double Ap = 0.;
      if (strstr (s, "ne") != NULL)
	Ap = 0.;
      else if (strstr (s, "te") != NULL)
	Ap = ap;
      else if (strstr (s, "ni") != NULL)
	Ap = 0.;
      else if (strstr (s, "ti") != NULL)
	Ap = ap;
      else if (strstr (s, "nb") != NULL)
	Ap = 0.;
      else if (strstr (s, "pb") != NULL)
	Ap = 0.;
      else if (strstr (s, "ptot") != NULL)
	Ap = ap;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	Ap = 0.;
      else if (strstr (s, "omegp") != NULL)
	Ap = 0.;
      else if (strstr (s, "omgvb") != NULL)
	Ap = 0.;
      else if (strstr (s, "omgpp") != NULL)
	Ap = 0.;
      else if (strstr (s, "omgeb") != NULL)
	Ap = 0.;
      else if (strstr (s, "er") != NULL)
	Ap = 0.;
      else if (strstr (s, "ommvb") != NULL)
	Ap = 0.;
      else if (strstr (s, "ommpp") != NULL)
	Ap = 0.;
      else if (strstr (s, "omevb") != NULL)
	Ap = 0.;
      else if (strstr (s, "omepp") != NULL)
	Ap = 0.;
      else if (strstr (s, "kpol") != NULL)
	Ap = 0.;
      else if (strstr (s, "omghb") != NULL)
	Ap = 0.;
      else if (strstr (s, "nz1") != NULL)
	Ap = 0.;
      else if (strstr (s, "vtor1") != NULL)
	Ap = 0.;
      else if (strstr (s, "vpol1") != NULL)
	Ap = 0.;

      if (strstr (s, "te") != NULL)
	field.ShiftScale (ine, Ap);
      else if (strstr (s, "ti") != NULL)
	field.ShiftScale (ini, Ap);
      else
	field.Shift (Ap);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));

      printf ("pFile field: %-34s:  shift = %11.4e\n", s, Ap);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 6 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType6 (double aw, Field& RBp, Field& R)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType6: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType6: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType6: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double Aw = 0.;
      if (strstr (s, "ne") != NULL)
	Aw = 0.;
      else if (strstr (s, "te") != NULL)
	Aw = 0.;
      else if (strstr (s, "ni") != NULL)
	Aw = 0.;
      else if (strstr (s, "ti") != NULL)
	Aw = 0.;
      else if (strstr (s, "nb") != NULL)
	Aw = 0.;
      else if (strstr (s, "pb") != NULL)
	Aw = 0.;
      else if (strstr (s, "ptot") != NULL)
	Aw = 0.;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	Aw = aw;
      else if (strstr (s, "omegp") != NULL)
	Aw = 0.;
      else if (strstr (s, "omgvb") != NULL)
	Aw = aw;
      else if (strstr (s, "omgpp") != NULL)
	Aw = 0.;
      else if (strstr (s, "omgeb") != NULL)
	Aw = aw;
      else if (strstr (s, "er") != NULL)
	Aw = aw;
      else if (strstr (s, "ommvb") != NULL)
	Aw = aw;
      else if (strstr (s, "ommpp") != NULL)
	Aw = 0.;
      else if (strstr (s, "omevb") != NULL)
	Aw = aw;
      else if (strstr (s, "omepp") != NULL)
	Aw = 0.;
      else if (strstr (s, "kpol") != NULL)
	Aw = 0.;
      else if (strstr (s, "omghb") != NULL)
	Aw = 0.;
      else if (strstr (s, "nz1") != NULL)
	Aw = 0.;
      else if (strstr (s, "vtor1") != NULL)
	Aw = aw;
      else if (strstr (s, "vpol1") != NULL)
	Aw = 0.;

      if (strstr (s, "er") != NULL)
	field.ShiftScale (RBp, Aw);
      else if (strstr (s, "vtor1") != NULL)
	field.ShiftScale (R, Aw);
      else
	field.Shift (Aw);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));

      printf ("pFile field: %-34s:  shift = %11.4e\n", s, Aw);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 7 rescaling of pFile 
// #############################################
void Rescale::pFileRescaleType7 (double a1)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];
  Field  field;
  
  // Check for existence of input pfile
  FILE* file = OpenFiler ((char*) "Inputs/pFile");
  if (file == NULL) 
    {
      printf ("RESCALE::pFileRescaleType7: Error opening input pFile\n");
      exit (1);
    }

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/pFile");
 
  do
    {
      // Read field
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("RESCALE::pFileRescaleType7: Error reading input pFile\n");
	  exit (1);
	}
      
      field.resize (n);
      for (int i = 0; i < n; i++)
	{
	  if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
	    {
	      printf ("RESCALE::pFileRescaleType7: Error reading %s\n", s);
	      exit (1);
	    }
	  else
	    {
	      field.PushData (i, x, y, dydx);
	    }
	}

      // Rescale field
      double A1 = 1.;
      if (strstr (s, "ne") != NULL)
	A1 = a1;
      else if (strstr (s, "te") != NULL)
	A1  = a1;
      else if (strstr (s, "ni") != NULL)
	A1 = a1;
      else if (strstr (s, "ti") != NULL)
	A1  = a1;
      else if (strstr (s, "nb") != NULL)
	A1 = a1;
      else if (strstr (s, "pb") != NULL)
	A1 = a1*a1;
      else if (strstr (s, "ptot") != NULL)
	A1 = a1*a1;
      else if (strstr (s, "omeg") != NULL && strstr (s, "omegp") == NULL)
	A1 = 1.;
      else if (strstr (s, "omegp") != NULL)
	A1 = 1.;
      else if (strstr (s, "omgvb") != NULL)
	A1 = 1.;
      else if (strstr (s, "omgpp") != NULL)
	A1 = 1.;
      else if (strstr (s, "omgeb") != NULL)
	A1 = 1.;
      else if (strstr (s, "er") != NULL)
	A1 = a1;
      else if (strstr (s, "ommvb") != NULL)
	A1 = 1.;
      else if (strstr (s, "ommpp") != NULL)
	A1 = 1.;
      else if (strstr (s, "omevb") != NULL)
	A1 = 1.;
      else if (strstr (s, "omepp") != NULL)
	A1 = 1.;
      else if (strstr (s, "kpol") != NULL)
	A1 = 1.;
      else if (strstr (s, "omghb") != NULL)
	A1 = a1;
      else if (strstr (s, "nz1") != NULL)
	A1 = a1;
      else if (strstr (s, "vtor1") != NULL)
	A1 = 1.;
      else if (strstr (s, "vpol1") != NULL)
	A1 = a1;

      field.Rescale (A1);

      // Output field to output pFile
      fprintf (file1, "%d %s\n", n, s);
      for (int i = 0; i < n; i++)
	fprintf (file1, "%e %e %e\n", field.GetX(i), field.GetY(i), field.GetdYdX(i));
      

      printf ("pFile field: %-34s:  rescale factor = %11.4e\n", s, A1);
    }
  while (1);

  fclose (file); fclose (file1);
}

// #############################################
// Function to perform Type 4 rescaling of cFile
// #############################################
void Rescale::cFileRescaleType4 (double ac)
{
  int    n;
  double x, y1, y2, y3, y4, dydx;
  Field  Chip, Chie, Chin, Chii;
    
  // Check for existence of cFile
  FILE* file = OpenFiler ((char*) "Inputs/cFile");
  if (file == NULL) 
    {
      printf ("RESCALE::cFileRescaleType4: Error opening cFile\n");
      exit (1);
    }

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("RESCALE::cFileRescaleType4: Error reading cFile (1)\n");
      exit (1);
    }

  Chip.resize (n);
  Chie.resize (n);
  Chin.resize (n);
  Chii.resize (n);

  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("RESCALE::cFileRescaleType4: Error reading cFile (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip.PushData (i, x, y1, dydx);
	  Chie.PushData (i, x, y2, dydx);
	  Chin.PushData (i, x, y3, dydx);
	  Chii.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // Open output pfile
  FILE* file1 = OpenFilew ((char*) "Outputs/cFile");
  fprintf (file1, "%d\n", n);
  for (int i = 0; i < n; i++)
    fprintf (file, "%19.6e %19.6e %19.6e %19.6e %19.6e\n", Chip.GetX (i), ac * Chip.GetY (i), ac * Chie.GetY (i), ac * Chin.GetY (i), ac * Chii.GetY (i));

  fclose (file1);

  printf ("Type 4 cFile rescaling:  scale factor = %11.4e\n", ac);
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Rescale::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("RESCALE::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for reading
// ##########################################
FILE* Rescale::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("RESCALE::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ############################################
// Function to open existing file for appending
// ############################################
FILE* Rescale::OpenFilea (char* filename)
{
  FILE* file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("RESCALE::OpenFilea: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to call operating system
// #################################
void Rescale::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("RESCALE: Operating system call error executing %s\n", command);
      exit (1);
    }
}

