// cFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Neoclassical:: cFileInterp               (vector<string> cFileName,   vector<double> cFileTime,   int cFileNumber, double time)
// void Neoclassical:: cFileInterpolateLinear    (char* cFile1, double time1, char* cFile2, double time2, char* cFile,     double time)
// void Neoclassical:: cFileInterpolateQuadratic (char* cFile1, double time1, char* cFile2, double time2, char* cFile,     double time)
// void Neoclassical:: cFileInterpolateCubic     (char* cFile1, double time1, char* cFile2, double time2, char* cFile3,    double time3, char* cFile, double time)
// void Neoclassical:: cFileInterpolateQuartic   (char* cFile1, double time1, char* cFile2, double time2, char* cFile3,    double time3,
//				                  char* cFile4, double time4, char* cFile,  double time)

#include "Neoclassical.h"

// ###############################
// Functions to interpolate cFiles
// ###############################
void Neoclassical::cFileInterp (vector<string> cFileName, vector<double> cFileTime, int cFileNumber, double time)
{
  if (cFileNumber < 1)
    {
      printf ("NEOCLASSICAL::cFileInterp - cFileNumber must be greater than zero\n");
      exit (1);
    }
  else if (cFileNumber == 1)
    {
      char* cFile = "Inputs/cFile";
      char* file1 = (char*) cFileName[0].c_str();

      cFileInterpolateLinear (file1, cFileTime[0], cFile, time);
    }
  else if (cFileNumber == 2)
    {
      char* cFile = "Inputs/cFile";
      char* file1 = (char*) cFileName[0].c_str();
      char* file2 = (char*) cFileName[1].c_str();

      cFileInterpolateQuadratic (file1, cFileTime[0], file2, cFileTime[1], cFile, time);
    }
  else if (CATS)
    {
      int index;

      if (time < cFileTime[0])
	index = 0;
      else if (time >= cFileTime[cFileNumber-1])
	index = cFileNumber - 2;
      else
	{
	  for (int i = 0; i < cFileNumber - 1; i++)
	    if (time >= cFileTime[i] && time < cFileTime[i+1])
	      index = i;
	}
      
      char* cFile = "Inputs/cFile";
      char* file1 = (char*) cFileName[index  ].c_str();
      char* file2 = (char*) cFileName[index+1].c_str();

      cFileInterpolateQuadratic (file1, cFileTime[index], file2, cFileTime[index+1], cFile, time);
    }
  else if (cFileNumber == 3)
    {
      char* cFile = "Inputs/cFile";
      char* file1 = (char*) cFileName[0].c_str();
      char* file2 = (char*) cFileName[1].c_str();
      char* file3 = (char*) cFileName[2].c_str();

      cFileInterpolateCubic (file1, cFileTime[0], file2, cFileTime[1], file3, cFileTime[2], cFile, time);
    }
  else if (cFileNumber == 4)
    {
      char* cFile = "Inputs/cFile";
      char* file1 = (char*) cFileName[0].c_str();
      char* file2 = (char*) cFileName[1].c_str();
      char* file3 = (char*) cFileName[2].c_str();
      char* file4 = (char*) cFileName[3].c_str();

      cFileInterpolateQuartic (file1, cFileTime[0], file2, cFileTime[1], file3, cFileTime[2],
			       file4, cFileTime[3], cFile, time);
    }
  else
    {
      int index, cntrl;

      if (time < cFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= cFileTime[cFileNumber-1])
	{
	  index = cFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < cFileNumber-1; i++)
	    if (time >= cFileTime[i] && time < cFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == cFileNumber-2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* cFile = "Inputs/cFile";
	  char* file1 = (char*) cFileName[index-1].c_str();
	  char* file2 = (char*) cFileName[index  ].c_str();
	  char* file3 = (char*) cFileName[index+1].c_str();
	  char* file4 = (char*) cFileName[index+2].c_str();
	  
	  cFileInterpolateQuartic (file1, cFileTime[index-1], file2, cFileTime[index], file3, cFileTime[index+1],
				   file4, cFileTime[index+2], cFile, time);
	}
      else if (cntrl == 2)
	{
	  char* cFile = "Inputs/cFile";
	  char* file1 = (char*) cFileName[index  ].c_str();
	  char* file2 = (char*) cFileName[index+1].c_str();
	  char* file3 = (char*) cFileName[index+2].c_str();
	  
	  cFileInterpolateCubic (file1, cFileTime[index], file2, cFileTime[index+1], file3, cFileTime[index+2], cFile, time);
	}
      else if (cntrl == 3)
	{
	  char* cFile = "Inputs/cFile";
	  char* file1 = (char*) cFileName[index-1].c_str();
	  char* file2 = (char*) cFileName[index  ].c_str();
	  char* file3 = (char*) cFileName[index+1].c_str();
	  
	  cFileInterpolateCubic (file1, cFileTime[index-1], file2, cFileTime[index], file3, cFileTime[index+1], cFile, time);
	}
    }
}

void Neoclassical::cFileInterpolateLinear (char* cFile1, double time1, char* cFile, double time)
{
  int    n;
  double x, y1, y2, y3, y4, dydx;

  // ################
  // Read first cFile
  // ################
  Field Chip_1, Chie_1, Chin_1, Chii_1;        

  FILE* file = OpenFiler (cFile1);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateLinear: Error reading cFile_1 (1)\n");
      exit (1);
    }

  Chip_1.resize (n);
  Chie_1.resize (n);
  Chin_1.resize (n);
  Chii_1.resize (n);

  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateLinear: Error reading cFile_1 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_1.PushData (i, x, y1, dydx);
	  Chie_1.PushData (i, x, y2, dydx);
	  Chin_1.PushData (i, x, y3, dydx);
	  Chii_1.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // ######################
  // Interpolate cFile data
  // ######################
  double weight1 = 1.;
  
  Field Chip, Chie, Chin, Chii;
  FieldInterpolateLinear (Chip_1, Chip, weight1);
  FieldInterpolateLinear (Chie_1, Chie, weight1);
  FieldInterpolateLinear (Chin_1, Chin, weight1);
  FieldInterpolateLinear (Chii_1, Chii, weight1);
  
  // ########################
  // Write interpolated cFile
  // ########################
  file = OpenFilew (cFile);

  fprintf (file, "%d\n", Chip.GetN ());
  for (int i = 0; i < Chip.GetN (); i++)
    {
      Chip.PullData (i, x, y1, dydx);
      Chie.PullData (i, x, y2, dydx);
      Chin.PullData (i, x, y3, dydx);
      Chii.PullData (i, x, y4, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", x, y1, y2, y3, y4);
    }
    
  fclose (file);

  printf ("cFile Interpolation:\n");
  printf ("%s %11.4e\n", cFile1, weight1);
}

void Neoclassical::cFileInterpolateQuadratic (char* cFile1, double time1, char* cFile2, double time2, char* cFile, double time)
{
  int    n;
  double x, y1, y2, y3, y4, dydx;

  // ################
  // Read first cFile
  // ################
  Field Chip_1, Chie_1, Chin_1, Chii_1;         

  FILE* file = OpenFiler (cFile1);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateQuadratic: Error reading cFile_1 (1)\n");
      exit (1);
    }

  Chip_1.resize (n);
  Chie_1.resize (n);
  Chin_1.resize (n);
  Chii_1.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateQuadratic: Error reading cFile_1 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_1.PushData (i, x, y1, dydx);
	  Chie_1.PushData (i, x, y2, dydx);
	  Chin_1.PushData (i, x, y3, dydx);
	  Chii_1.PushData (i, x, y4, dydx);

	}
    }
  
  fclose (file);

  // #################
  // Read second cFile
  // #################
  Field Chip_2, Chie_2, Chin_2, Chii_2;        

  file = OpenFiler (cFile2);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateQuadratic: Error reading cFile_2 (1)\n");
      exit (1);
    }

  Chip_2.resize (n);
  Chie_2.resize (n);
  Chin_2.resize (n);
  Chii_2.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateQuadratic: Error reading cFile_2 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_2.PushData (i, x, y1, dydx);
	  Chie_2.PushData (i, x, y2, dydx);
	  Chin_2.PushData (i, x, y3, dydx);
	  Chii_2.PushData (i, x, y4, dydx);

	}
    }
  
  fclose (file);

  // ######################
  // Interpolate cFile data
  // ######################
  double weight1 = (time - time2) /(time1 - time2);
  double weight2 = (time - time1) /(time2 - time1);
  
  Field Chip, Chie, Chin, Chii;
  FieldInterpolateQuadratic (Chip_1, Chip_2, Chip, weight1, weight2);
  FieldInterpolateQuadratic (Chie_1, Chie_2, Chie, weight1, weight2);
  FieldInterpolateQuadratic (Chin_1, Chin_2, Chin, weight1, weight2);
  FieldInterpolateQuadratic (Chii_1, Chii_2, Chii, weight1, weight2);
  
  // ########################
  // Write interpolated cFile
  // ########################
  file = OpenFilew (cFile);

  fprintf (file, "%d\n", Chip.GetN ());
  for (int i = 0; i < Chip.GetN (); i++)
    {
      Chip.PullData (i, x, y1, dydx);
      Chie.PullData (i, x, y2, dydx);
      Chin.PullData (i, x, y3, dydx);
      Chii.PullData (i, x, y4, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", x, y1, y2, y3, y4);
    }
    
  fclose (file);

  printf ("cFile Interpolation:\n");
  printf ("%s %11.4e\n", cFile1, weight1);
  printf ("%s %11.4e\n", cFile2, weight2);
}

void Neoclassical::cFileInterpolateCubic (char* cFile1, double time1, char* cFile2, double time2, char* cFile3, double time3, char* cFile, double time)
{
  int    n;
  double x, y1, y2, y3, y4, dydx;

  // ################
  // Read first cFile
  // ################
  Field Chip_1, Chie_1, Chin_1, Chii_1;  

  FILE* file = OpenFiler (cFile1);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_1 (1)\n");
      exit (1);
    }

  Chip_1.resize (n);
  Chie_1.resize (n);
  Chin_1.resize (n);
  Chii_1.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_1 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_1.PushData (i, x, y1, dydx);
	  Chie_1.PushData (i, x, y2, dydx);
	  Chin_1.PushData (i, x, y3, dydx);
	  Chii_1.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // #################
  // Read second cFile
  // #################
  Field Chip_2, Chie_2, Chin_2, Chii_2;          

  file = OpenFiler (cFile2);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_2 (1)\n");
      exit (1);
    }

  Chip_2.resize (n);
  Chie_2.resize (n);
  Chin_2.resize (n);
  Chii_2.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_2 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_2.PushData (i, x, y1, dydx);
	  Chie_2.PushData (i, x, y2, dydx);
	  Chin_2.PushData (i, x, y3, dydx);
	  Chii_2.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // ################
  // Read third cFile
  // ################
  Field Chip_3, Chie_3, Chin_3, Chii_3;

  file = OpenFiler (cFile3);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_3 (1)\n");
      exit (1);
    }

  Chip_3.resize (n);
  Chie_3.resize (n);
  Chin_3.resize (n);
  Chii_3.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_3 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_3.PushData (i, x, y1, dydx);
	  Chie_3.PushData (i, x, y2, dydx);
	  Chin_3.PushData (i, x, y3, dydx);
	  Chii_3.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // #######################
  // Interpolate cFile data
  // #######################
  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);

  Field Chip, Chie, Chin, Chii;
  FieldInterpolateCubic (Chip_1, Chip_2, Chip_3, Chip, weight1, weight2, weight3);
  FieldInterpolateCubic (Chie_1, Chie_2, Chie_3, Chie, weight1, weight2, weight3);
  FieldInterpolateCubic (Chin_1, Chin_2, Chin_3, Chin, weight1, weight2, weight3);
  FieldInterpolateCubic (Chii_1, Chii_2, Chii_3, Chii, weight1, weight2, weight3);
 
  // ########################
  // Write interpolated cFile
  // ########################
  file = OpenFilew (cFile);

  fprintf (file, "%d\n", Chip.GetN ());
  for (int i = 0; i < Chip.GetN (); i++)
    {
      Chip.PullData (i, x, y1, dydx);
      Chie.PullData (i, x, y2, dydx);
      Chin.PullData (i, x, y3, dydx);
      Chii.PullData (i, x, y4, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", x, y1, y2, y3, y4);
    }
    
  fclose (file);

  printf ("cFile Interpolation:\n");
  printf ("%s %11.4e\n", cFile1, weight1);
  printf ("%s %11.4e\n", cFile2, weight2);
  printf ("%s %11.4e\n", cFile3, weight3);
}

void Neoclassical::cFileInterpolateQuartic (char* cFile1, double time1, char* cFile2, double time2, char* cFile3, double time3, char* cFile4,
					    double time4, char* cFile, double time)
{
  int    n;
  double x, y1, y2, y3, y4, dydx;

  // ################
  // Read first cFile
  // ################
  Field Chip_1, Chie_1, Chin_1, Chii_1; 

  FILE* file = OpenFiler (cFile1);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_1 (1)\n");
      exit (1);
    }

  Chip_1.resize (n);
  Chie_1.resize (n);
  Chin_1.resize (n);
  Chii_1.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_1 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_1.PushData (i, x, y1, dydx);
	  Chie_1.PushData (i, x, y2, dydx);
	  Chin_1.PushData (i, x, y3, dydx);
	  Chii_1.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // #################
  // Read second cFile
  // #################
  Field Chip_2, Chie_2, Chin_2, Chii_2;       

  file = OpenFiler (cFile2);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_2 (1)\n");
      exit (1);
    }

  Chip_2.resize (n);
  Chie_2.resize (n);
  Chin_2.resize (n);
  Chii_2.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_2 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_2.PushData (i, x, y1, dydx);
	  Chie_2.PushData (i, x, y2, dydx);
	  Chin_2.PushData (i, x, y3, dydx);
	  Chii_2.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // ################
  // Read third cFile
  // ################
  Field Chip_3, Chie_3, Chin_3, Chii_3;        

  file = OpenFiler (cFile3);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_3 (1)\n");
      exit (1);
    }

  Chip_3.resize (n);
  Chie_3.resize (n);
  Chin_3.resize (n);
  Chii_3.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_3 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_3.PushData (i, x, y1, dydx);
	  Chie_3.PushData (i, x, y2, dydx);
	  Chin_3.PushData (i, x, y3, dydx);
	  Chii_3.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // #################
  // Read fourth cFile
  // #################
  Field Chip_4, Chie_4, Chin_4, Chii_4;        

  file = OpenFiler (cFile4);

  // Read data from cFile
  if (fscanf (file, "%d", &n) != 1)
    {
      printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_4 (1)\n");
      exit (1);
    }

  Chip_4.resize (n);
  Chie_4.resize (n);
  Chin_4.resize (n);
  Chii_4.resize (n);
  
  for (int i = 0; i < n; i++)
    {
      if (fscanf (file, "%lf %lf %lf %lf %lf", &x, &y1, &y2, &y3, &y4) != 5)
	{
	  printf ("NEOCLASSICAL::cFileInterpolateCubic: Error reading cFile_4 (2)\n");
	  exit (1);
	}
      else
	{
	  dydx = 0.;
	  Chip_4.PushData (i, x, y1, dydx);
	  Chie_4.PushData (i, x, y2, dydx);
	  Chin_4.PushData (i, x, y3, dydx);
	  Chii_4.PushData (i, x, y4, dydx);
	}
    }
  
  fclose (file);

  // ######################
  // Interpolate cFile data
  // ######################
  double weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4);
  double weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4);
  double weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4);
  double weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3);

  Field Chip, Chie, Chin, Chii;
  FieldInterpolateQuartic (Chip_1, Chip_2, Chip_3, Chip_3, Chip, weight1, weight2, weight3, weight4);
  FieldInterpolateQuartic (Chie_1, Chie_2, Chie_3, Chie_3, Chie, weight1, weight2, weight3, weight4);
  FieldInterpolateQuartic (Chin_1, Chin_2, Chin_3, Chin_3, Chin, weight1, weight2, weight3, weight4);
  FieldInterpolateQuartic (Chii_1, Chii_2, Chii_3, Chii_3, Chii, weight1, weight2, weight3, weight4);

  // ########################
  // Write interpolated cFile
  // ########################
  file = OpenFilew (cFile);

  fprintf (file, "%d\n", Chip.GetN ());
  for (int i = 0; i < Chip.GetN (); i++)
    {
      Chip.PullData (i, x, y1, dydx);
      Chie.PullData (i, x, y2, dydx);
      Chin.PullData (i, x, y3, dydx);
      Chii.PullData (i, x, y4, dydx);
      fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e\n", x, y1, y2, y3, y4);
    }
    
  fclose (file);
  
  printf ("cFile Interpolation:\n");
  printf ("%s %11.4e\n", cFile1, weight1);
  printf ("%s %11.4e\n", cFile2, weight2);
  printf ("%s %11.4e\n", cFile3, weight3);
  printf ("%s %11.4e\n", cFile4, weight4);
}

