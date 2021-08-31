// pFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Neoclassical:: pFileInterp               (vector<string> pFileName,   vector<double> pFileTime,   int pFileNumber,  double time)
// void Neoclassical:: pFileInterpolateLinear    (char* pFile1, double time1, char* pFile, double time)
// void Neoclassical:: pFileInterpolateQuadratic (char* pFile1, double time1, char* pFile2, double time2, char* pFile,      double time)
// void Neoclassical:: pFileInterpolateCubic     (char* pFile1, double time1, char* pFile2, double time2, char* pFile3,     double time3, char* pFile, double time)
// void Neoclassical:: pFileInterpolateQuartic   (char* pFile1, double time1, char* pFile2, double time2, char* pFile3,     double time3,
//				                  char* pFile4, double time4, char* pFile,  double time)
// void Neoclassical:: FieldInterpolateLinear    (Field& Field1, Field& Field,  double weight1)
// void Neoclassical:: FieldInterpolateQuadratic (Field& Field1, Field& Field2, Field& Field, double weight1, double weight2)
// void Neoclassical:: FieldInterpolateCubic     (Field& Field1, Field& Field2, Field& Field3, Field& Field, double weight1, double weight2, double weight3)
// void Neoclassical:: FieldInterpolateQuartic   (Field& Field1, Field& Field2, Field& Field3, Field& Field4, Field& Field,
//					          double weight1, double weight2, double weight3, double weight4)

#include "Neoclassical.h"

// ###############################
// Functions to interpolate pFiles
// ###############################
void Neoclassical::pFileInterp (vector<string> pFileName, vector<double> pFileTime, int pFileNumber, double time)
{
  if (pFileNumber < 1)
    {
      printf ("NEOCLASSICAL::pFileInterp - pFileNumber must be greater than zero\n");
      exit (1);
    }
  else if (pFileNumber == 1)
    {
      char* pFile = "Inputs/pFile";
      char* file1 = (char*) pFileName[0].c_str();
 
      pFileInterpolateLinear (file1, pFileTime[0], pFile, time);
    }
  else if (pFileNumber == 2)
    {
      char* pFile = "Inputs/pFile";
      char* file1 = (char*) pFileName[0].c_str();
      char* file2 = (char*) pFileName[1].c_str();

      pFileInterpolateQuadratic (file1, pFileTime[0], file2, pFileTime[1], pFile, time);
    }
  else if (pFileNumber == 3)
    {
      char* pFile = "Inputs/pFile";
      char* file1 = (char*) pFileName[0].c_str();
      char* file2 = (char*) pFileName[1].c_str();
      char* file3 = (char*) pFileName[2].c_str();

      pFileInterpolateCubic (file1, pFileTime[0], file2, pFileTime[1], file3, pFileTime[2], pFile, time);
    }
  else if (pFileNumber == 4)
    {
      char* pFile = "Inputs/pFile";
      char* file1 = (char*) pFileName[0].c_str();
      char* file2 = (char*) pFileName[1].c_str();
      char* file3 = (char*) pFileName[2].c_str();
      char* file4 = (char*) pFileName[3].c_str();

      pFileInterpolateQuartic (file1, pFileTime[0], file2, pFileTime[1], file3, pFileTime[2],
			       file4, pFileTime[3], pFile, time);
    }
  else
    {
      int index, cntrl;

      if (time < pFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= pFileTime[pFileNumber-1])
	{
	  index = pFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < pFileNumber-1; i++)
	    if (time >= pFileTime[i] && time < pFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == pFileNumber-2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* pFile = "Inputs/pFile";
	  char* file1 = (char*) pFileName[index-1].c_str();
	  char* file2 = (char*) pFileName[index  ].c_str();
	  char* file3 = (char*) pFileName[index+1].c_str();
	  char* file4 = (char*) pFileName[index+2].c_str();
	  
	  pFileInterpolateQuartic (file1, pFileTime[index-1], file2, pFileTime[index], file3, pFileTime[index+1],
				   file4, pFileTime[index+2], pFile, time);
	}
      else if (cntrl == 2)
	{
	  char* pFile = "Inputs/pFile";
	  char* file1 = (char*) pFileName[index  ].c_str();
	  char* file2 = (char*) pFileName[index+1].c_str();
	  char* file3 = (char*) pFileName[index+2].c_str();
	  
	  pFileInterpolateCubic (file1, pFileTime[index], file2, pFileTime[index+1], file3, pFileTime[index+2], pFile, time);
	}
      else if (cntrl == 3)
	{
	  char* pFile = "Inputs/pFile";
	  char* file1 = (char*) pFileName[index-1].c_str();
	  char* file2 = (char*) pFileName[index  ].c_str();
	  char* file3 = (char*) pFileName[index+1].c_str();
	  
	  pFileInterpolateCubic (file1, pFileTime[index-1], file2, pFileTime[index], file3, pFileTime[index+1], pFile, time);
	}
    }
}

void Neoclassical::pFileInterpolateLinear (char* pFile1, double time1, char* pFile, double time)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];

  // ################
  // Read first pFile
  // ################
 
  Field ne_1;          Field Te_1;        Field ni_1;        Field Ti_1;  Field wp_1;
  Field nb_1;          Field wE_1;        Field nI_1;        Field NZA_1; Field wt_1;
  int   ne_flag_1 = 0; int Te_flag_1 = 0; int ni_flag_1 = 0; int Ti_flag_1 = 0;
  int   nb_flag_1 = 0; int wE_flag_1 = 0; int nI_flag_1 = 0; int NZ_flag_1 = 0;
  int   wt_flag_1 = 0; int wp_flag_1 = 0;

  FILE* file = OpenFiler (pFile1);
 
  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading pFile_1\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_1 = 1;
	  printf ("Reading ne    from pFile_1 - n = %4d:\n", n);
	  ne_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_1 = 1;
	  printf ("Reading Te    from pFile_1 - n = %4d:\n", n);
	  Te_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_1 = 1;
	  printf ("Reading ni    from pFile_1 - n = %4d:\n", n);
	  ni_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_1 = 1;
	  printf ("Reading Ti    from pFile_1 - n = %4d:\n", n);
	  Ti_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_1 = 1;
	  printf ("Reading nb    from pFile_1 - n = %4d:\n", n);
	  nb_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_1 = 1;
	  printf ("Reading omgeb from pFile_1 - n = %4d:\n", n);
	  wE_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_1 = 1;
	  printf ("Reading  omeg from pFile_1 - n = %4d:\n", n);
	  wt_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_1.PushData (i, x, y, dydx);
		}
	    }
	}
       else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_1 = 1;
	  printf ("Reading omegp from pFile_1 - n = %4d:\n", n);
	  wp_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_1 = 1;
	  printf ("Reading nz1   from pFile_1 - n = %4d:\n", n);
	  nI_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_1 = 1;
	  printf ("Reading NZA   from pFile_1 - n = %4d:\n", n);
	  NZA_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateLinear: Error reading unused field in pFile_1\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_1 * Te_flag_1 * ni_flag_1 * Ti_flag_1 * nb_flag_1 * wE_flag_1 * nI_flag_1 * NZ_flag_1 * wt_flag_1 * wp_flag_1 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateLinear: Missing field in pFile_1\n");
      exit (1);
    }

  // ########################
  // Interpolate profile data
  // ########################
  double weight1 = 1.;
  
  Field ne;
  FieldInterpolateLinear (ne_1,  ne,  weight1);
  Field Te;
  FieldInterpolateLinear (Te_1,  Te,  weight1);
  Field ni;
  FieldInterpolateLinear (ni_1,  ni,  weight1);
  Field Ti;
  FieldInterpolateLinear (Ti_1,  Ti,  weight1);
  Field nb;
  FieldInterpolateLinear (nb_1,  nb,  weight1);
  Field wE;
  FieldInterpolateLinear (wE_1,  wE,  weight1);
  Field wt;
  FieldInterpolateLinear (wt_1,  wt,  weight1);
  Field wp;
  FieldInterpolateLinear (wp_1,  wp,  weight1);
  Field nI;
  FieldInterpolateLinear (nI_1,  nI,  weight1);
  Field NZA;
  FieldInterpolateLinear (NZA_1, NZA, weight1);

  // ########################
  // Write interpolated pFile
  // ########################
  file = OpenFilew (pFile);

  fprintf (file, "%3d %s %s %s", ne.GetN (), "PsiN", " ne(10^20/m^3)", " dne/dpsiN\n");
  for (int i = 0; i < ne.GetN (); i++)
    {
      ne.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Te.GetN (), "PsiN", " te(KeV)", " dte/dpsiN\n");
  for (int i = 0; i < Te.GetN (); i++)
    {
      Te.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ni.GetN (), "PsiN", " ni(10^20/m^3)", " dni/dpsiN\n");
  for (int i = 0; i < ni.GetN (); i++)
    {
      ni.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Ti.GetN (), "PsiN", " ti(keV)", " dti/dpsiN\n");
  for (int i = 0; i < Ti.GetN (); i++)
    {
      Ti.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nb.GetN (), "PsiN", " nb(10^20/m^3)", " dnb/dpsiN\n");
  for (int i = 0; i < nb.GetN (); i++)
    {
      nb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wE.GetN (), "PsiN", " omgeb(kRad/s)", " domgeb/dpsiN\n");
  for (int i = 0; i < wE.GetN (); i++)
    {
      wE.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wt.GetN (), "PsiN", " omeg(kRad/s)", " domeg/dpsiN\n");
  for (int i = 0; i < wt.GetN (); i++)
    {
      wt.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wp.GetN (), "PsiN", " omegp(kRad/s)", " domegp/dpsiN\n");
  for (int i = 0; i < wp.GetN (); i++)
    {
      wp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nI.GetN (), "PsiN", " nz1(10^20/m^3)", " dnz1/dpsiN\n");
  for (int i = 0; i < nI.GetN (); i++)
    {
      nI.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", NZA.GetN (), "N", " Z", " A of ION SPECIES\n");
  for (int i = 0; i < NZA.GetN (); i++)
    {
      NZA.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
    
  fclose (file);

  printf ("pFile Interpolation:\n");
  printf ("%s %11.4e\n", pFile1, weight1);
}

void Neoclassical::pFileInterpolateQuadratic (char* pFile1, double time1, char* pFile2, double time2, char* pFile, double time)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];

  // ################
  // Read first pFile
  // ################
 
  Field ne_1;          Field Te_1;        Field ni_1;        Field Ti_1;  Field wp_1;
  Field nb_1;          Field wE_1;        Field nI_1;        Field NZA_1; Field wt_1;
  int   ne_flag_1 = 0; int Te_flag_1 = 0; int ni_flag_1 = 0; int Ti_flag_1 = 0;
  int   nb_flag_1 = 0; int wE_flag_1 = 0; int nI_flag_1 = 0; int NZ_flag_1 = 0;
  int   wt_flag_1 = 0; int wp_flag_1 = 0;
  
  FILE* file = OpenFiler (pFile1);
 
  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading pFile_1\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_1 = 1;
	  printf ("Reading ne    from pFile_1 - n = %4d:\n", n);
	  ne_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_1 = 1;
	  printf ("Reading Te    from pFile_1 - n = %4d:\n", n);
	  Te_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_1 = 1;
	  printf ("Reading ni    from pFile_1 - n = %4d:\n", n);
	  ni_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_1 = 1;
	  printf ("Reading Ti    from pFile_1 - n = %4d:\n", n);
	  Ti_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_1 = 1;
	  printf ("Reading nb    from pFile_1 - n = %4d:\n", n);
	  nb_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_1 = 1;
	  printf ("Reading omgeb from pFile_1 - n = %4d:\n", n);
	  wE_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_1 = 1;
	  printf ("Reading  omeg from pFile_1 - n = %4d:\n", n);
	  wt_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_1 = 1;
	  printf ("Reading omegp from pFile_1 - n = %4d:\n", n);
	  wp_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_1 = 1;
	  printf ("Reading nz1   from pFile_1 - n = %4d:\n", n);
	  nI_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_1 = 1;
	  printf ("Reading NZA   from pFile_1 - n = %4d:\n", n);
	  NZA_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading unused field in pFile_1\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_1 * Te_flag_1 * ni_flag_1 * Ti_flag_1 * nb_flag_1 * wE_flag_1 * nI_flag_1 * NZ_flag_1 * wt_flag_1 * wp_flag_1 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Missing field in pFile_1\n");
      exit (1);
    }

  // #################
  // Read second pFile
  // #################
 
  Field ne_2;          Field Te_2;        Field ni_2;        Field Ti_2;  Field wp_2;
  Field nb_2;          Field wE_2;        Field nI_2;        Field NZA_2; Field wt_2;
  int   ne_flag_2 = 0; int Te_flag_2 = 0; int ni_flag_2 = 0; int Ti_flag_2 = 0;
  int   nb_flag_2 = 0; int wE_flag_2 = 0; int nI_flag_2 = 0; int NZ_flag_2 = 0;
  int   wt_flag_2 = 0; int wp_flag_2 = 0;
  
  file = OpenFiler (pFile2);

  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading pFile_2\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_2 = 1;
	  printf ("Reading ne    from pFile_2 - n = %4d:\n", n);
	  ne_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_2 = 1;
	  printf ("Reading Te    from pFile_2 - n = %4d:\n", n);
	  Te_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_2 = 1;
	  printf ("Reading ni    from pFile_2 - n = %4d:\n", n);
	  ni_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_2 = 1;
	  printf ("Reading Ti    from pFile_2 - n = %4d:\n", n);
	  Ti_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_2 = 1;
	  printf ("Reading nb    from pFile_2 - n = %4d:\n", n);
	  nb_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_2 = 1;
	  printf ("Reading omgeb from pFile_2 - n = %4d:\n", n);
	  wE_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_2 = 1;
	  printf ("Reading  omeg from pFile_2 - n = %4d:\n", n);
	  wt_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_2.PushData (i, x, y, dydx);
		}
	    }
	}
       else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_2 = 1;
	  printf ("Reading omegp from pFile_2 - n = %4d:\n", n);
	  wp_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_2 = 1;
	  printf ("Reading nz1   from pFile_2 - n = %4d:\n", n);
	  nI_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_2 = 1;
	  printf ("Reading NZA   from pFile_2 - n = %4d:\n", n);
	  NZA_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Error reading unused field in pFile_2\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_2 * Te_flag_2 * ni_flag_2 * Ti_flag_2 * nb_flag_2 * wE_flag_2 * nI_flag_2 * NZ_flag_2 * wt_flag_2 * wp_flag_2 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateQuadratic: Missing field in pFile_2\n");
      exit (1);
    }

  // ########################
  // Interpolate profile data
  // ########################
  double weight1 = (time - time2) /(time1 - time2);
  double weight2 = (time - time1) /(time2 - time1);
  
  Field ne;
  FieldInterpolateQuadratic (ne_1,  ne_2,  ne,  weight1, weight2);
  Field Te;
  FieldInterpolateQuadratic (Te_1,  Te_2,  Te,  weight1, weight2);
  Field ni;
  FieldInterpolateQuadratic (ni_1,  ni_2,  ni,  weight1, weight2);
  Field Ti;
  FieldInterpolateQuadratic (Ti_1,  Ti_2,  Ti,  weight1, weight2);
  Field nb;
  FieldInterpolateQuadratic (nb_1,  nb_2,  nb,  weight1, weight2);
  Field wE;
  FieldInterpolateQuadratic (wE_1,  wE_2,  wE,  weight1, weight2);
  Field wt;
  FieldInterpolateQuadratic (wt_1,  wt_2,  wt,  weight1, weight2);
  Field wp;
  FieldInterpolateQuadratic (wp_1,  wp_2,  wp,  weight1, weight2);
  Field nI;
  FieldInterpolateQuadratic (nI_1,  nI_2,  nI,  weight1, weight2);
  Field NZA;
  FieldInterpolateQuadratic (NZA_1, NZA_2, NZA, weight1, weight2);

  // ########################
  // Write interpolated pFile
  // ########################
  file = OpenFilew (pFile);

  fprintf (file, "%3d %s %s %s", ne.GetN (), "PsiN", " ne(10^20/m^3)", " dne/dpsiN\n");
  for (int i = 0; i < ne.GetN (); i++)
    {
      ne.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Te.GetN (), "PsiN", " te(KeV)", " dte/dpsiN\n");
  for (int i = 0; i < Te.GetN (); i++)
    {
      Te.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ni.GetN (), "PsiN", " ni(10^20/m^3)", " dni/dpsiN\n");
  for (int i = 0; i < ni.GetN (); i++)
    {
      ni.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Ti.GetN (), "PsiN", " ti(keV)", " dti/dpsiN\n");
  for (int i = 0; i < Ti.GetN (); i++)
    {
      Ti.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nb.GetN (), "PsiN", " nb(10^20/m^3)", " dnb/dpsiN\n");
  for (int i = 0; i < nb.GetN (); i++)
    {
      nb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wE.GetN (), "PsiN", " omgeb(kRad/s)", " domgeb/dpsiN\n");
  for (int i = 0; i < wE.GetN (); i++)
    {
      wE.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wt.GetN (), "PsiN", " omeg(kRad/s)", " domeg/dpsiN\n");
  for (int i = 0; i < wt.GetN (); i++)
    {
      wt.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wp.GetN (), "PsiN", " omegp(kRad/s)", " domegp/dpsiN\n");
  for (int i = 0; i < wp.GetN (); i++)
    {
      wp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nI.GetN (), "PsiN", " nz1(10^20/m^3)", " dnz1/dpsiN\n");
  for (int i = 0; i < nI.GetN (); i++)
    {
      nI.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", NZA.GetN (), "N", " Z", " A of ION SPECIES\n");
  for (int i = 0; i < NZA.GetN (); i++)
    {
      NZA.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
    
  fclose (file);

  printf ("pFile Interpolation:\n");
  printf ("%s %11.4e\n", pFile1, weight1);
  printf ("%s %11.4e\n", pFile2, weight2);
}

void Neoclassical::pFileInterpolateCubic (char* pFile1, double time1, char* pFile2, double time2, char* pFile3, double time3, char* pFile, double time)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];

  // ################
  // Read first pFile
  // ################
 
  Field ne_1;          Field Te_1;        Field ni_1;        Field Ti_1;  Field wp_1;
  Field nb_1;          Field wE_1;        Field nI_1;        Field NZA_1; Field wt_1;
  int   ne_flag_1 = 0; int Te_flag_1 = 0; int ni_flag_1 = 0; int Ti_flag_1 = 0;
  int   nb_flag_1 = 0; int wE_flag_1 = 0; int nI_flag_1 = 0; int NZ_flag_1 = 0;
  int   wt_flag_1 = 0; int wp_flag_1 = 0;
  
  FILE* file = OpenFiler (pFile1);
 
  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading pFile_1\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_1 = 1;
	  printf ("Reading ne    from pFile_1 - n = %4d:\n", n);
	  ne_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_1 = 1;
	  printf ("Reading Te    from pFile_1 - n = %4d:\n", n);
	  Te_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_1 = 1;
	  printf ("Reading ni    from pFile_1 - n = %4d:\n", n);
	  ni_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_1 = 1;
	  printf ("Reading Ti    from pFile_1 - n = %4d:\n", n);
	  Ti_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_1 = 1;
	  printf ("Reading nb    from pFile_1 - n = %4d:\n", n);
	  nb_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_1 = 1;
	  printf ("Reading omgeb from pFile_1 - n = %4d:\n", n);
	  wE_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_1 = 1;
	  printf ("Reading omegp from pFile_1 - n = %4d:\n", n);
	  wp_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_1 = 1;
	  printf ("Reading  omeg from pFile_1 - n = %4d:\n", n);
	  wt_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_1 = 1;
	  printf ("Reading nz1   from pFile_1 - n = %4d:\n", n);
	  nI_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_1 = 1;
	  printf ("Reading NZA   from pFile_1 - n = %4d:\n", n);
	  NZA_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading unused field in pFile_1\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_1 * Te_flag_1 * ni_flag_1 * Ti_flag_1 * nb_flag_1 * wE_flag_1 * nI_flag_1 * NZ_flag_1 * wt_flag_1 * wp_flag_1 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateCubic: Missing field in pFile_1\n");
      exit (1);
    }

  // #################
  // Read second pFile
  // #################
 
  Field ne_2;          Field Te_2;        Field ni_2;        Field Ti_2;  Field wp_2;
  Field nb_2;          Field wE_2;        Field nI_2;        Field NZA_2; Field wt_2;
  int   ne_flag_2 = 0; int Te_flag_2 = 0; int ni_flag_2 = 0; int Ti_flag_2 = 0;
  int   nb_flag_2 = 0; int wE_flag_2 = 0; int nI_flag_2 = 0; int NZ_flag_2 = 0;
  int   wt_flag_2 = 0; int wp_flag_2 = 0;
  
  file = OpenFiler (pFile2);

  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading pFile_2\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_2 = 1;
	  printf ("Reading ne    from pFile_2 - n = %4d:\n", n);
	  ne_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_2 = 1;
	  printf ("Reading Te    from pFile_2 - n = %4d:\n", n);
	  Te_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_2 = 1;
	  printf ("Reading ni    from pFile_2 - n = %4d:\n", n);
	  ni_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_2 = 1;
	  printf ("Reading Ti    from pFile_2 - n = %4d:\n", n);
	  Ti_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_2 = 1;
	  printf ("Reading nb    from pFile_2 - n = %4d:\n", n);
	  nb_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_2 = 1;
	  printf ("Reading omgeb from pFile_2 - n = %4d:\n", n);
	  wE_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_2 = 1;
	  printf ("Reading  omeg from pFile_2 - n = %4d:\n", n);
	  wt_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_2 = 1;
	  printf ("Reading omegp from pFile_2 - n = %4d:\n", n);
	  wp_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_2 = 1;
	  printf ("Reading nz1   from pFile_2 - n = %4d:\n", n);
	  nI_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_2 = 1;
	  printf ("Reading NZA   from pFile_2 - n = %4d:\n", n);
	  NZA_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading unused field in pFile_2\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_2 * Te_flag_2 * ni_flag_2 * Ti_flag_2 * nb_flag_2 * wE_flag_2 * nI_flag_2 * NZ_flag_2 * wt_flag_2 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateCubic: Missing field in pFile_2\n");
      exit (1);
    }

  // ################
  // Read third pFile
  // ################
 
  Field ne_3;          Field Te_3;        Field ni_3;        Field Ti_3;  Field wp_3;
  Field nb_3;          Field wE_3;        Field nI_3;        Field NZA_3; Field wt_3;
  int   ne_flag_3 = 0; int Te_flag_3 = 0; int ni_flag_3 = 0; int Ti_flag_3 = 0;
  int   nb_flag_3 = 0; int wE_flag_3 = 0; int nI_flag_3 = 0; int NZ_flag_3 = 0;
  int   wt_flag_3 = 0; int wp_flag_3 = 0;
  
  file = OpenFiler (pFile3);

  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading pFile_3\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_3 = 1;
	  printf ("Reading ne    from pFile_3 - n = %4d:\n", n);
	  ne_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_3 = 1;
	  printf ("Reading Te    from pFile_3 - n = %4d:\n", n);
	  Te_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_3 = 1;
	  printf ("Reading ni    from pFile_3 - n = %4d:\n", n);
	  ni_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_3 = 1;
	  printf ("Reading Ti    from pFile_3 - n = %4d:\n", n);
	  Ti_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_3 = 1;
	  printf ("Reading nb    from pFile_3 - n = %4d:\n", n);
	  nb_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_3 = 1;
	  printf ("Reading omgeb from pFile_3 - n = %4d:\n", n);
	  wE_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_3 = 1;
	  printf ("Reading  omeg from pFile_3 - n = %4d:\n", n);
	  wt_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_3 = 1;
	  printf ("Reading omegp from pFile_3 - n = %4d:\n", n);
	  wp_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_3 = 1;
	  printf ("Reading nz1   from pFile_3 - n = %4d:\n", n);
	  nI_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_3 = 1;
	  printf ("Reading NZA   from pFile_3 - n = %4d:\n", n);
	  NZA_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateCubic: Error reading unused field in pFile_3\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_3 * Te_flag_3 * ni_flag_3 * Ti_flag_3 * nb_flag_3 * wE_flag_3 * nI_flag_3 * NZ_flag_3 * wt_flag_3 * wp_flag_3 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateCubic: Missing field in pFile_3\n");
      exit (1);
    }
  
  fclose (file);

  // ########################
  // Interpolate profile data
  // ########################
  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);

  Field ne;
  FieldInterpolateCubic (ne_1,  ne_2,  ne_3,  ne,  weight1, weight2, weight3);
  Field Te;
  FieldInterpolateCubic (Te_1,  Te_2,  Te_3,  Te,  weight1, weight2, weight3);
  Field ni;
  FieldInterpolateCubic (ni_1,  ni_2,  ni_3,  ni,  weight1, weight2, weight3);
  Field Ti;
  FieldInterpolateCubic (Ti_1,  Ti_2,  Ti_3,  Ti,  weight1, weight2, weight3);
  Field nb;
  FieldInterpolateCubic (nb_1,  nb_2,  nb_3,  nb,  weight1, weight2, weight3);
  Field wE;
  FieldInterpolateCubic (wE_1,  wE_2,  wE_3,  wE,  weight1, weight2, weight3);
  Field wt;
  FieldInterpolateCubic (wt_1,  wt_2,  wt_3,  wt,  weight1, weight2, weight3);
  Field wp;
  FieldInterpolateCubic (wp_1,  wp_2,  wp_3,  wp,  weight1, weight2, weight3);
  Field nI;
  FieldInterpolateCubic (nI_1,  nI_2,  nI_3,  nI,  weight1, weight2, weight3);
  Field NZA;
  FieldInterpolateCubic (NZA_1, NZA_2, NZA_3, NZA, weight1, weight2, weight3);

  // ########################
  // Write interpolated pFile
  // ########################
  file = OpenFilew (pFile);

  fprintf (file, "%3d %s %s %s", ne.GetN (), "PsiN", " ne(10^20/m^3)", " dne/dpsiN\n");
  for (int i = 0; i < ne.GetN (); i++)
    {
      ne.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Te.GetN (), "PsiN", " te(KeV)", " dte/dpsiN\n");
  for (int i = 0; i < Te.GetN (); i++)
    {
      Te.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ni.GetN (), "PsiN", " ni(10^20/m^3)", " dni/dpsiN\n");
  for (int i = 0; i < ni.GetN (); i++)
    {
      ni.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Ti.GetN (), "PsiN", " ti(keV)", " dti/dpsiN\n");
  for (int i = 0; i < Ti.GetN (); i++)
    {
      Ti.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nb.GetN (), "PsiN", " nb(10^20/m^3)", " dnb/dpsiN\n");
  for (int i = 0; i < nb.GetN (); i++)
    {
      nb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wE.GetN (), "PsiN", " omgeb(kRad/s)", " domgeb/dpsiN\n");
  for (int i = 0; i < wE.GetN (); i++)
    {
      wE.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wt.GetN (), "PsiN", " omeg(kRad/s)", " domeg/dpsiN\n");
  for (int i = 0; i < wt.GetN (); i++)
    {
      wt.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wp.GetN (), "PsiN", " omegp(kRad/s)", " domegp/dpsiN\n");
  for (int i = 0; i < wp.GetN (); i++)
    {
      wp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nI.GetN (), "PsiN", " nz1(10^20/m^3)", " dnz1/dpsiN\n");
  for (int i = 0; i < nI.GetN (); i++)
    {
      nI.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", NZA.GetN (), "N", " Z", " A of ION SPECIES\n");
  for (int i = 0; i < NZA.GetN (); i++)
    {
      NZA.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
    
  fclose (file);

  printf ("pFile Interpolation:\n");
  printf ("%s %11.4e\n", pFile1, weight1);
  printf ("%s %11.4e\n", pFile2, weight2);
  printf ("%s %11.4e\n", pFile3, weight3);
}

void Neoclassical::pFileInterpolateQuartic (char* pFile1, double time1, char* pFile2, double time2, char* pFile3, double time3, char* pFile4,
					    double time4, char* pFile, double time)
{
  int    n;
  double x, y, dydx;
  char   s[MAXFILENAMELENGTH];

  // ################
  // Read first pFile
  // ################
 
  Field ne_1;          Field Te_1;        Field ni_1;        Field Ti_1;  Field wp_1;
  Field nb_1;          Field wE_1;        Field nI_1;        Field NZA_1; Field wt_1;
  int   ne_flag_1 = 0; int Te_flag_1 = 0; int ni_flag_1 = 0; int Ti_flag_1 = 0;
  int   nb_flag_1 = 0; int wE_flag_1 = 0; int nI_flag_1 = 0; int NZ_flag_1 = 0;
  int   wt_flag_1 = 0; int wp_flag_1 = 0;
  
  FILE* file = OpenFiler (pFile1);
 
  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading pFile_1\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_1 = 1;
	  printf ("Reading ne    from pFile_1 - n = %4d:\n", n);
	  ne_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_1 = 1;
	  printf ("Reading Te    from pFile_1 - n = %4d:\n", n);
	  Te_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_1 = 1;
	  printf ("Reading ni    from pFile_1 - n = %4d:\n", n);
	  ni_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_1 = 1;
	  printf ("Reading Ti    from pFile_1 - n = %4d:\n", n);
	  Ti_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_1 = 1;
	  printf ("Reading nb    from pFile_1 - n = %4d:\n", n);
	  nb_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_1 = 1;
	  printf ("Reading omgeb from pFile_1 - n = %4d:\n", n);
	  wE_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_1 = 1;
	  printf ("Reading  omeg from pFile_1 - n = %4d:\n", n);
	  wt_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omep(") != NULL)
	{
	  // Read omep field (assumed units krad/s)
	  wp_flag_1 = 1;
	  printf ("Reading omegp from pFile_1 - n = %4d:\n", n);
	  wp_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_1 = 1;
	  printf ("Reading nz1   from pFile_1 - n = %4d:\n", n);
	  nI_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_1 = 1;
	  printf ("Reading NZA   from pFile_1 - n = %4d:\n", n);
	  NZA_1.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_1.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading unused field in pFile_1\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_1 * Te_flag_1 * ni_flag_1 * Ti_flag_1 * nb_flag_1 * wE_flag_1 * nI_flag_1 * NZ_flag_1 * wt_flag_1 * wp_flag_1 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateQuartic: Missing field in pFile_1\n");
      exit (1);
    }

  // #################
  // Read second pFile
  // #################
 
  Field ne_2;          Field Te_2;        Field ni_2;        Field Ti_2;  Field wp_2;
  Field nb_2;          Field wE_2;        Field nI_2;        Field NZA_2; Field wt_2;
  int   ne_flag_2 = 0; int Te_flag_2 = 0; int ni_flag_2 = 0; int Ti_flag_2 = 0;
  int   nb_flag_2 = 0; int wE_flag_2 = 0; int nI_flag_2 = 0; int NZ_flag_2 = 0;
  int   wt_flag_2 = 0; int wp_flag_2 = 0;
    
  file = OpenFiler (pFile2);

  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading pFile_2\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_2 = 1;
	  printf ("Reading ne    from pFile_2 - n = %4d:\n", n);
	  ne_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_2 = 1;
	  printf ("Reading Te    from pFile_2 - n = %4d:\n", n);
	  Te_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_2 = 1;
	  printf ("Reading ni    from pFile_2 - n = %4d:\n", n);
	  ni_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_2 = 1;
	  printf ("Reading Ti    from pFile_2 - n = %4d:\n", n);
	  Ti_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_2 = 1;
	  printf ("Reading nb    from pFile_2 - n = %4d:\n", n);
	  nb_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_2 = 1;
	  printf ("Reading omgeb from pFile_2 - n = %4d:\n", n);
	  wE_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_2 = 1;
	  printf ("Reading  omeg from pFile_2 - n = %4d:\n", n);
	  wt_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_2 = 1;
	  printf ("Reading omegp from pFile_2 - n = %4d:\n", n);
	  wp_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_2 = 1;
	  printf ("Reading nz1   from pFile_2 - n = %4d:\n", n);
	  nI_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_2 = 1;
	  printf ("Reading NZA   from pFile_2 - n = %4d:\n", n);
	  NZA_2.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_2.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading unused field in pFile_2\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_2 * Te_flag_2 * ni_flag_2 * Ti_flag_2 * nb_flag_2 * wE_flag_2 * nI_flag_2 * NZ_flag_2 * wt_flag_2 * wp_flag_2 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateQuartic: Missing field in pFile_2\n");
      exit (1);
    }

  // ################
  // Read third pFile
  // ################
 
  Field ne_3;          Field Te_3;        Field ni_3;        Field Ti_3;  Field wp_3;
  Field nb_3;          Field wE_3;        Field nI_3;        Field NZA_3; Field wt_3;
  int   ne_flag_3 = 0; int Te_flag_3 = 0; int ni_flag_3 = 0; int Ti_flag_3 = 0;
  int   nb_flag_3 = 0; int wE_flag_3 = 0; int nI_flag_3 = 0; int NZ_flag_3 = 0;
  int   wt_flag_3 = 0; int wp_flag_3 = 0;
  
  file = OpenFiler (pFile3);

  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading pFile_3\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_3 = 1;
	  printf ("Reading ne    from pFile_3 - n = %4d:\n", n);
	  ne_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_3 = 1;
	  printf ("Reading Te    from pFile_3 - n = %4d:\n", n);
	  Te_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_3 = 1;
	  printf ("Reading ni    from pFile_3 - n = %4d:\n", n);
	  ni_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_3 = 1;
	  printf ("Reading Ti    from pFile_3 - n = %4d:\n", n);
	  Ti_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_3 = 1;
	  printf ("Reading nb    from pFile_3 - n = %4d:\n", n);
	  nb_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_3 = 1;
	  printf ("Reading omgeb from pFile_3 - n = %4d:\n", n);
	  wE_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_3 = 1;
	  printf ("Reading  omeg from pFile_3 - n = %4d:\n", n);
	  wt_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omegp(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_3 = 1;
	  printf ("Reading omegp from pFile_3 - n = %4d:\n", n);
	  wp_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_3 = 1;
	  printf ("Reading nz1   from pFile_3 - n = %4d:\n", n);
	  nI_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_3 = 1;
	  printf ("Reading NZA   from pFile_3 - n = %4d:\n", n);
	  NZA_3.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_3.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading unused field in pFile_3\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_3 * Te_flag_3 * ni_flag_3 * Ti_flag_3 * nb_flag_3 * wE_flag_3 * nI_flag_3 * NZ_flag_3 * wt_flag_3 * wp_flag_3 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateQuartic: Missing field in pFile_3\n");
      exit (1);
    }
  
  fclose (file);

  // #################
  // Read fourth pFile
  // #################
 
  Field ne_4;          Field Te_4;        Field ni_4;        Field Ti_4;  Field wp_4;
  Field nb_4;          Field wE_4;        Field nI_4;        Field NZA_4; Field wt_4;
  int   ne_flag_4 = 0; int Te_flag_4 = 0; int ni_flag_4 = 0; int Ti_flag_4 = 0;
  int   nb_flag_4 = 0; int wE_flag_4 = 0; int nI_flag_4 = 0; int NZ_flag_4 = 0;
  int   wt_flag_4 = 0; int wp_flag_4 = 0;
  
  file = OpenFiler (pFile4);

  do
    {
      int retval = fscanf (file, "%d %[^\n]s", &n, &s);

      if (retval == EOF)
	break;
      
      if (retval != 2)
	{
	  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading pFile_4\n");
	  exit (1);
	}

      if (strstr (s, "ne") != NULL)
	{
	  // Read ne field (assumed units - 10^20 m^-3)
	  ne_flag_4 = 1;
	  printf ("Reading ne    from pFile_4 - n = %4d:\n", n);
	  ne_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ne\n");
		  exit (1);
		}
	      else
		{
		  ne_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "te") != NULL)
	{
	  // Read Te field (assumed units - keV)
	  Te_flag_4 = 1;
	  printf ("Reading Te    from pFile_4 - n = %4d:\n", n);
	  Te_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Te\n");
		  exit (1);
		}
	      else
		{
		  Te_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ni") != NULL)
	{
	  // Read ni field (assumed units - 10^20 m^-3)
	  ni_flag_4 = 1;
	  printf ("Reading ni    from pFile_4 - n = %4d:\n", n);
	  ni_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading ni\n");
		  exit (1);
		}
	      else
		{
		  ni_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ti") != NULL)
	{
	  // Read Ti field (assumed units - keV)
	  Ti_flag_4 = 1;
	  printf ("Reading Ti    from pFile_4 - n = %4d:\n", n);
	  Ti_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading Ti\n");
		  exit (1);
		}
	      else
		{
		  Ti_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nb") != NULL)
	{
	  // Read nb field (assumed units - 10^20 m^-3)
	  nb_flag_4 = 1;
	  printf ("Reading nb    from pFile_4 - n = %4d:\n", n);
	  nb_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nb\n");
		  exit (1);
		}
	      else
		{
		  nb_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omgeb") != NULL)
	{
	  // Read omgeb field (assumed units krad/s)
	  wE_flag_4 = 1;
	  printf ("Reading omgeb from pFile_4 - n = %4d:\n", n);
	  wE_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omgeb\n");
		  exit (1);
		}
	      else
		{
		  wE_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omeg field (assumed units krad/s)
	  wt_flag_4 = 1;
	  printf ("Reading  omeg from pFile_4 - n = %4d:\n", n);
	  wt_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omeg\n");
		  exit (1);
		}
	      else
		{
		  wt_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "omeg(") != NULL)
	{
	  // Read omegp field (assumed units krad/s)
	  wp_flag_4 = 1;
	  printf ("Reading omegp from pFile_4 - n = %4d:\n", n);
	  wp_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading omegp\n");
		  exit (1);
		}
	      else
		{
		  wp_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "nz1") != NULL)
	{
	  // Read nz1 field (assumed units 10^20 m^-3)
	  nI_flag_4 = 1;
	  printf ("Reading nz1   from pFile_4 - n = %4d:\n", n);
	  nI_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading nI1\n");
		  exit (1);
		}
	      else
		{
		  nI_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else if (strstr (s, "ION") != NULL)
	{
	  // Read NZA field
	  NZ_flag_4 = 1;
	  printf ("Reading NZA   from pFile_4 - n = %4d:\n", n);
	  NZA_4.resize (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading NZA\n");
		  exit (1);
		}
	      else
		{
		  NZA_4.PushData (i, x, y, dydx);
		}
	    }
	}
      else
	{
	  // Read unused fields
	  Field Unused (n);
	  for (int i = 0; i < n; i++)
	    {
	      if (fscanf (file, "%lf %lf %lf", &x, &y, &dydx) != 3)
		{
		  printf ("NEOCLASSICAL::pFileInterpolateQuartic: Error reading unused field in pFile_4\n");
		  exit (1);
		}
	      else
		{
		  Unused.PushData (i, x, y, dydx);
		}
	    }
	}
    }
  while (1);
  
  fclose (file);

  if (ne_flag_4 * Te_flag_4 * ni_flag_4 * Ti_flag_4 * nb_flag_4 * wE_flag_4 * nI_flag_4 * NZ_flag_4 * wt_flag_4 * wp_flag_4 == 0)
    {
      printf ("NEOCLASSICAL::pFileInterpolateQuartic: Missing field in pFile_4\n");
      exit (1);
    }
  
  fclose (file);

  // ########################
  // Interpolate profile data
  // ########################
  double weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4);
  double weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4);
  double weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4);
  double weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3);

  Field ne;
  FieldInterpolateQuartic (ne_1,  ne_2,  ne_3,  ne_4,  ne,  weight1, weight2, weight3, weight4);
  Field Te;
  FieldInterpolateQuartic (Te_1,  Te_2,  Te_3,  Te_4,  Te,  weight1, weight2, weight3, weight4);
  Field ni;
  FieldInterpolateQuartic (ni_1,  ni_2,  ni_3,  ni_4,  ni,  weight1, weight2, weight3, weight4);
  Field Ti;
  FieldInterpolateQuartic (Ti_1,  Ti_2,  Ti_3,  Ti_4,  Ti,  weight1, weight2, weight3, weight4);
  Field nb;
  FieldInterpolateQuartic (nb_1,  nb_2,  nb_3,  nb_4,  nb,  weight1, weight2, weight3, weight4);
  Field wE;
  FieldInterpolateQuartic (wE_1,  wE_2,  wE_3,  wE_4,  wE,  weight1, weight2, weight3, weight4);
  Field wt;
  FieldInterpolateQuartic (wt_1,  wt_2,  wt_3,  wt_4,  wt,  weight1, weight2, weight3, weight4);
  Field wp;
  FieldInterpolateQuartic (wp_1,  wp_2,  wp_3,  wp_4,  wp,  weight1, weight2, weight3, weight4);
  Field nI;
  FieldInterpolateQuartic (nI_1,  nI_2,  nI_3,  nI_4,  nI,  weight1, weight2, weight3, weight4);
  Field NZA;
  FieldInterpolateQuartic (NZA_1, NZA_2, NZA_3, NZA_4, NZA, weight1, weight2, weight3, weight4);

  // ########################
  // Write interpolated pFile
  // ########################
  file = OpenFilew (pFile);

  fprintf (file, "%3d %s %s %s", ne.GetN (), "PsiN", " ne(10^20/m^3)", " dne/dpsiN\n");
  for (int i = 0; i < ne.GetN (); i++)
    {
      ne.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Te.GetN (), "PsiN", " te(KeV)", " dte/dpsiN\n");
  for (int i = 0; i < Te.GetN (); i++)
    {
      Te.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", ni.GetN (), "PsiN", " ni(10^20/m^3)", " dni/dpsiN\n");
  for (int i = 0; i < ni.GetN (); i++)
    {
      ni.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", Ti.GetN (), "PsiN", " ti(keV)", " dti/dpsiN\n");
  for (int i = 0; i < Ti.GetN (); i++)
    {
      Ti.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nb.GetN (), "PsiN", " nb(10^20/m^3)", " dnb/dpsiN\n");
  for (int i = 0; i < nb.GetN (); i++)
    {
      nb.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wE.GetN (), "PsiN", " omgeb(kRad/s)", " domgeb/dpsiN\n");
  for (int i = 0; i < wE.GetN (); i++)
    {
      wE.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wt.GetN (), "PsiN", " omeg(kRad/s)", " domeg/dpsiN\n");
  for (int i = 0; i < wt.GetN (); i++)
    {
      wt.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", wp.GetN (), "PsiN", " omegp(kRad/s)", " domegp/dpsiN\n");
  for (int i = 0; i < wp.GetN (); i++)
    {
      wp.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", nI.GetN (), "PsiN", " nz1(10^20/m^3)", " dnz1/dpsiN\n");
  for (int i = 0; i < nI.GetN (); i++)
    {
      nI.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
  fprintf (file, "%3d %s %s %s", NZA.GetN (), "N", " Z", " A of ION SPECIES\n");
  for (int i = 0; i < NZA.GetN (); i++)
    {
      NZA.PullData (i, x, y, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e\n", x, y, dydx);
    }
    
  fclose (file);

  printf ("pFile Interpolation:\n");
  printf ("%s %11.4e\n", pFile1, weight1);
  printf ("%s %11.4e\n", pFile2, weight2);
  printf ("%s %11.4e\n", pFile3, weight3);
  printf ("%s %11.4e\n", pFile4, weight4);
}

// ###############################
// Functions to interpolate Fields
// ###############################
void Neoclassical::FieldInterpolateLinear (Field& Field1, Field& Field, double weight1)
{
  int n;
  int n1 = Field1.GetN ();

  n = n1;

  Field.resize (n);

  double x1, y1, dydx1;
  double x,  y,  dydx;
  for (int i = 0; i < n; i++)
    {
      Field1.PullData (i, x1, y1, dydx1);

      x    = x1;
      y    = y1;
      dydx = dydx1;

      Field.PushData (i, x, y, dydx);
    }
}

void Neoclassical::FieldInterpolateQuadratic (Field& Field1, Field& Field2, Field& Field, double weight1, double weight2)
{
  int n;
  int n1 = Field1.GetN ();
  int n2 = Field2.GetN ();

  if (n1 == n2)
    {
      n = n1;
    }
  else
    {
      printf ("NEOCLASSICAL::FieldInterpolateQuadratic: Error - Field size mismatch\n");
      exit (1);
    }

  Field.resize (n);

  double x1, y1, dydx1;
  double x2, y2, dydx2;
  double x,  y,  dydx;
  for (int i = 0; i < n; i++)
    {
      Field1.PullData (i, x1, y1, dydx1);
      Field2.PullData (i, x2, y2, dydx2);

      x    = weight1 * x1    + weight2 * x2;
      y    = weight1 * y1    + weight2 * y2;
      dydx = weight1 * dydx1 + weight2 * dydx2;

      Field.PushData (i, x, y, dydx);
    }
}

void Neoclassical::FieldInterpolateCubic (Field& Field1, Field& Field2, Field& Field3, Field& Field, double weight1, double weight2, double weight3)
{
  int n;
  int n1 = Field1.GetN ();
  int n2 = Field2.GetN ();
  int n3 = Field3.GetN ();

  if (n1 == n2 && n2 == n3)
    {
      n = n1;
    }
  else
    {
      printf ("NEOCLASSICAL::FieldInterpolateCubic: Error - Field size mismatch\n");
      exit (1);
    }

  Field.resize (n);

  double x1, y1, dydx1;
  double x2, y2, dydx2;
  double x3, y3, dydx3;
  double x,  y,  dydx;
  for (int i = 0; i < n; i++)
    {
      Field1.PullData (i, x1, y1, dydx1);
      Field2.PullData (i, x2, y2, dydx2);
      Field3.PullData (i, x3, y3, dydx3);

      x    = weight1 * x1    + weight2 * x2    + weight3 * x3;
      y    = weight1 * y1    + weight2 * y2    + weight3 * y3;
      dydx = weight1 * dydx1 + weight2 * dydx2 + weight3 * dydx3;

      Field.PushData (i, x, y, dydx);
    }
}

void Neoclassical::FieldInterpolateQuartic (Field& Field1, Field& Field2, Field& Field3, Field& Field4, Field& Field,
					    double weight1, double weight2, double weight3, double weight4)
{
  int n;
  int n1 = Field1.GetN ();
  int n2 = Field2.GetN ();
  int n3 = Field3.GetN ();
  int n4 = Field4.GetN ();

  if (n1 == n2 && n2 == n3 && n3 == n4)
    {
      n = n1;
    }
  else
    {
      printf ("NEOCLASSICAL::FieldInterpolateQuartic: Error - Field size mismatch\n");
      exit (1);
    }

  Field.resize (n);

  double x1, y1, dydx1;
  double x2, y2, dydx2;
  double x3, y3, dydx3;
  double x4, y4, dydx4;
  double x,  y,  dydx;
  for (int i = 0; i < n; i++)
    {
      Field1.PullData (i, x1, y1, dydx1);
      Field2.PullData (i, x2, y2, dydx2);
      Field3.PullData (i, x3, y3, dydx3);
      Field4.PullData (i, x4, y4, dydx4);

      x    = weight1 * x1    + weight2 * x2    + weight3 * x3    + weight4 * x4;
      y    = weight1 * y1    + weight2 * y2    + weight3 * y3    + weight4 * y4;
      dydx = weight1 * dydx1 + weight2 * dydx2 + weight3 * dydx3 + weight4 * dydx4;

      Field.PushData (i, x, y, dydx);
    }
}
