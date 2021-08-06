// fFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Neoclassical:: fFileInterp               (vector<string> fFileName,   vector<double> fFileTime,   int fFileNumber,  double time)
// void Neoclassical:: fFileInterpolateLinear    (char* fFile1, double time1, char* fFile,  double time)
// void Neoclassical:: fFileInterpolateQuadratic (char* fFile1, double time1, char* fFile2, double time2, char* fFile,      double time)
// void Neoclassical:: fFileInterpolateCubic     (char* fFile1, double time1, char* fFile2, double time2, char* fFile3,     double time3, char* fFile, double time)
// void Neoclassical:: fFileInterpolateQuartic   (char* fFile1, double time1, char* fFile2, double time2, char* fFile3,     double time3,
//				                  char* fFile4, double time4, char* fFile,  double time)

#include "Neoclassical.h"

// ###############################
// Functions to interpolate fFiles
// ###############################
void Neoclassical::fFileInterp (vector<string> fFileName, vector<double> fFileTime, int fFileNumber, double time)
{
  if (fFileNumber < 1)
    {
      printf ("NEOCLASSICAL::fFileInterp - fFileNumber must be greater than zero\n");
      exit (1);
    }
  else if (fFileNumber == 1)
    {
      char* fFile = "Inputs/fFile";
      char* file1 = (char*) fFileName[0].c_str();

      fFileInterpolateLinear (file1, fFileTime[0], fFile, time);
    }
  else if (fFileNumber == 2)
    {
      char* fFile = "Inputs/fFile";
      char* file1 = (char*) fFileName[0].c_str();
      char* file2 = (char*) fFileName[1].c_str();

      fFileInterpolateQuadratic (file1, fFileTime[0], file2, fFileTime[1], fFile, time);
    }
  else if (fFileNumber == 3)
    {
      char* fFile = "Inputs/fFile";
      char* file1 = (char*) fFileName[0].c_str();
      char* file2 = (char*) fFileName[1].c_str();
      char* file3 = (char*) fFileName[2].c_str();

      fFileInterpolateCubic (file1, fFileTime[0], file2, fFileTime[1], file3, fFileTime[2], fFile, time);
    }
  else if (fFileNumber == 4)
    {
      char* fFile = "Inputs/fFile";
      char* file1 = (char*) fFileName[0].c_str();
      char* file2 = (char*) fFileName[1].c_str();
      char* file3 = (char*) fFileName[2].c_str();
      char* file4 = (char*) fFileName[3].c_str();

      fFileInterpolateQuartic (file1, fFileTime[0], file2, fFileTime[1], file3, fFileTime[2],
			       file4, fFileTime[3], fFile, time);
    }
  else
    {
      int index, cntrl = 0;

      if (time < fFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= fFileTime[fFileNumber-1])
	{
	  index = fFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < fFileNumber-1; i++)
	    if (time >= fFileTime[i] && time < fFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == fFileNumber-2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* fFile = "Inputs/fFile";
	  char* file1 = (char*) fFileName[index-1].c_str();
	  char* file2 = (char*) fFileName[index  ].c_str();
	  char* file3 = (char*) fFileName[index+1].c_str();
	  char* file4 = (char*) fFileName[index+2].c_str();
	  
	  fFileInterpolateQuartic (file1, fFileTime[index-1], file2, fFileTime[index], file3, fFileTime[index+1],
				   file4, fFileTime[index+2], fFile, time);
	}
      else if (cntrl == 2)
	{
	  char* fFile = "Inputs/fFile";
	  char* file1 = (char*) fFileName[index  ].c_str();
	  char* file2 = (char*) fFileName[index+1].c_str();
	  char* file3 = (char*) fFileName[index+2].c_str();
	  
	  fFileInterpolateCubic (file1, fFileTime[index], file2, fFileTime[index+1], file3, fFileTime[index+2], fFile, time);
	}
      else if (cntrl == 3)
	{
	  char* fFile = "Inputs/fFile";
	  char* file1 = (char*) fFileName[index-1].c_str();
	  char* file2 = (char*) fFileName[index  ].c_str();
	  char* file3 = (char*) fFileName[index+1].c_str();
	  
	  fFileInterpolateCubic (file1, fFileTime[index-1], file2, fFileTime[index], file3, fFileTime[index+1], fFile, time);
	}
      else
	{
	  printf ("NEOCLASSICAL::fFileInterp - Error cntrl = %1d\n", cntrl);
	  exit (1);
	}
    }
}

void Neoclassical::fFileInterpolateLinear (char* fFile1, double time1, char* fFile, double time)
{
  int ini; double inr;

  double          r1_1, r2_1, r3_1, r4_1, r5_1, r6_1, r7_1, r8_1, r9_1, r10_1, r11_1, r12_1, r13_1;
  int             NPSI_1, NTOR_1, nres_1;
  Array<double,1> v1_1, v2_1, v3_1;
  Array<int, 1>   mres_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1, u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1;
  Array<double,2> Freal_1, Fimag_1, Ereal_1, Eimag_1;
  
  double          r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, r10_0, r11_0, r12_0, r13_0;
  int             NPSI_0, NTOR_0, nres_0;
  Array<double,1> v1_0, v2_0, v3_0;
  Array<int, 1>   mres_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0, u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0;
  Array<double,2> Freal_0, Fimag_0, Ereal_0, Eimag_0;

  // ................
  // Read first fFile
  // ................
  FILE* file = OpenFiler (fFile1);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_1, &r2_1, &r3_1, &r4_1, &r5_1, &r6_1, &r7_1, &r8_1, &r9_1, &NPSI_1, &NTOR_1, &nres_1, &r10_1, &r11_1, &r12_1, &r13_1) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateLinear: Error reading fFile_1 (1)\n");
      exit (1);
    }

  v1_1.resize (NPSI_1); v2_1.resize (NPSI_1); v3_1.resize (NPSI_1);

  for (int j = 0; j < NPSI_1; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_1(j), &v2_1(j), &v3_1(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateLinear: Error reading fFile_1 (2)\n");
	  exit (1);
	}
    }

  mres_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize  (nres_1); u3_1.resize  (nres_1);
  u4_1.resize   (nres_1); u5_1.resize  (nres_1); u6_1.resize  (nres_1);
  u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1);
  u10_1.resize  (nres_1); u11_1.resize (nres_1); u12_1.resize (nres_1);
  u13_1.resize  (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize (nres_1);

  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &u1_1(j), &u2_1(j), &u3_1(j), &u4_1(j), &u5_1(j), &u6_1(j), &u7_1(j), &u8_1(j), &u9_1(j),
		  &u10_1(j), &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j), &u16_1(j), &u17_1(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateLinear: Error reading fFile_1 (3)\n");
	  exit (1);
	}
    }

  Freal_1.resize (nres_1, nres_1); Fimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_1(j, k), &Fimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateLinear: Error reading fFile_1 (4)\n");
	    exit (1);
	  }
      }
  
  Ereal_1.resize (nres_1, nres_1); Eimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_1(j, k), &Eimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateLinear: Error reading fFile_1 (5)\n");
	    exit (1);
	  }
      }
 
  // ......................
  // Interpolate fFile data
  // ......................

  NPSI_0 = NPSI_1;
  NTOR_0 = NTOR_1;
  nres_0 = nres_1;

  v1_0.resize     (NPSI_0); v2_0.resize  (NPSI_0); v3_0.resize  (NPSI_0);

  mres_0.resize   (nres_0);
  u1_0.resize     (nres_0); u2_0.resize  (nres_0); u3_0.resize  (nres_0);
  u4_0.resize     (nres_0); u5_0.resize  (nres_0); u6_0.resize  (nres_0);
  u7_0.resize     (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0);
  u10_0.resize    (nres_0); u11_0.resize (nres_0); u12_0.resize (nres_0);
  u13_0.resize    (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize    (nres_0); u17_0.resize (nres_0);

  Freal_0.resize  (nres_0, nres_0); Fimag_0.resize  (nres_0, nres_0);
  Ereal_0.resize  (nres_0, nres_0); Eimag_0.resize  (nres_0, nres_0);

  for (int i = 0; i < nres_0; i++)
    mres_0(i) = mres_1(i);
  
  double weight1 = 1.;
 
  r1_0  = r1_1;
  r2_0  = r2_1;
  r3_0  = r3_1; 
  r4_0  = r4_1; 
  r5_0  = r5_1; 
  r6_0  = r6_1; 
  r7_0  = r7_1; 
  r8_0  = r8_1; 
  r9_0  = r9_1;
  r10_0 = r10_1;
  r11_0 = r11_1;
  r12_0 = r12_1;
  r13_0 = r13_1;

  for (int j = 0; j < NPSI_0; j++)
    {
      v1_0(j) = v1_1(j);
      v2_0(j) = v2_1(j);
      v3_0(j) = v3_1(j);
    }

  for (int j = 0; j < nres_0; j++)
    {
      u1_0  (j) = u1_1  (j);
      u2_0  (j) = u2_1  (j);
      u3_0  (j) = u3_1  (j);
      u4_0  (j) = u4_1  (j);
      u5_0  (j) = u5_1  (j);
      u6_0  (j) = u6_1  (j);
      u7_0  (j) = u7_1  (j);
      u8_0  (j) = u8_1  (j);
      u9_0  (j) = u9_1  (j);
      u10_0 (j) = u10_1 (j);
      u11_0 (j) = u11_1 (j);
      u12_0 (j) = u12_1 (j);
      u13_0 (j) = u13_1 (j);
      u14_0 (j) = u14_1 (j);
      u15_0 (j) = u15_1 (j);
      u16_0 (j) = u16_1 (j);
      u17_0 (j) = u17_1 (j);
    }

  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      {
	Freal_0(j, k) = Freal_1(j, k);
	Fimag_0(j, k) = Fimag_1(j, k);
	Ereal_0(j, k) = Ereal_1(j, k);
	Eimag_0(j, k) = Eimag_1(j, k);
      }
    
  // ........................
  // Write interpolated fFile
  // ........................
  file = OpenFilew (fFile);
  
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e %16.9e %16.9e %16.9e\n",
	   r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, NPSI_0, NTOR_0, nres_0, r10_0, r11_0, r12_0, r13_0);
  for (int j = 0; j < NPSI_0; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     v1_0(j), v2_0(j), v3_0(j));
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j), u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Freal_0(j, k), Fimag_0(j, k));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Ereal_0(j, k), Eimag_0(j, k));
  
  fclose (file);
  
  printf ("fFile Interpolation:\n");
  printf ("%s %11.4e\n", fFile1, weight1);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "fFile Interpolation:\n");
  fprintf (monitor, "%s %11.4e\n", fFile1, weight1);
  fclose  (monitor);
}

void Neoclassical::fFileInterpolateQuadratic (char* fFile1, double time1, char* fFile2, double time2, char* fFile, double time)
{
  int ini; double inr;

  double          r1_1, r2_1, r3_1, r4_1, r5_1, r6_1, r7_1, r8_1, r9_1, r10_1, r11_1, r12_1, r13_1;
  int             NPSI_1, NTOR_1, nres_1;
  Array<double,1> v1_1, v2_1, v3_1;
  Array<int, 1>   mres_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1, u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1;
  Array<double,2> Freal_1, Fimag_1, Ereal_1, Eimag_1;
  
  double          r1_2, r2_2, r3_2, r4_2, r5_2, r6_2, r7_2, r8_2, r9_2, r10_2, r11_2, r12_2, r13_2;
  int             NPSI_2, NTOR_2, nres_2;
  Array<double,1> v1_2, v2_2, v3_2;
  Array<int, 1>   mres_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2, u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2;
  Array<double,2> Freal_2, Fimag_2, Ereal_2, Eimag_2;
 
  double          r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, r10_0, r11_0, r12_0, r13_0;
  int             NPSI_0, NTOR_0, nres_0;
  Array<double,1> v1_0, v2_0, v3_0;
  Array<int, 1>   mres_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0, u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0;
  Array<double,2> Freal_0, Fimag_0, Ereal_0, Eimag_0;

  // ................
  // Read first fFile
  // ................
  FILE* file = OpenFiler (fFile1);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_1, &r2_1, &r3_1, &r4_1, &r5_1, &r6_1, &r7_1, &r8_1, &r9_1, &NPSI_1, &NTOR_1, &nres_1, &r10_1, &r11_1, &r12_2, &r13_1) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_1 (1)\n");
      exit (1);
    }

  v1_1.resize (NPSI_1); v2_1.resize (NPSI_1); v3_1.resize (NPSI_1);

  for (int j = 0; j < NPSI_1; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_1(j), &v2_1(j), &v3_1(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_1 (2)\n");
	  exit (1);
	}
    }

  mres_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize  (nres_1); u3_1.resize  (nres_1);
  u4_1.resize   (nres_1); u5_1.resize  (nres_1); u6_1.resize  (nres_1);
  u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1);
  u10_1.resize  (nres_1); u11_1.resize (nres_1); u12_1.resize (nres_1);
  u13_1.resize  (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize (nres_1);
 
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &u1_1(j), &u2_1(j), &u3_1(j), &u4_1(j), &u5_1(j), &u6_1(j), &u7_1(j), &u8_1(j), &u9_1(j),
		  &u10_1(j), &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j), &u16_1(j), &u17_1(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_1 (3)\n");
	  exit (1);
	}
    }

  Freal_1.resize (nres_1, nres_1); Fimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_1(j, k), &Fimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_1 (4)\n");
	    exit (1);
	  }
      }
  
  Ereal_1.resize (nres_1, nres_1); Eimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_1(j, k), &Eimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_1 (5)\n");
	    exit (1);
	  }
      }
 
  fclose (file);

  // .................
  // Read second fFile
  // .................
  file = OpenFiler (fFile2);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_2, &r2_2, &r3_2, &r4_2, &r5_2, &r6_2, &r7_2, &r8_2, &r9_2, &NPSI_2, &NTOR_2, &nres_2, &r10_2, &r11_2, &r12_2, &r13_2) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_2 (1)\n");
      exit (1);
    }

  v1_2.resize (NPSI_2); v2_2.resize (NPSI_2); v3_2.resize (NPSI_2);

  for (int j = 0; j < NPSI_2; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_2(j), &v2_2(j), &v3_2(j)) != 3)
	{
	  printf ("'NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_2 (2)\n");
	  exit (1);
	}
    }

  mres_2.resize (nres_2);
  u1_2.resize   (nres_2); u2_2.resize  (nres_2); u3_2.resize  (nres_2);
  u4_2.resize   (nres_2); u5_2.resize  (nres_2); u6_2.resize  (nres_2);
  u7_2.resize   (nres_2); u8_2.resize  (nres_2); u9_2.resize  (nres_2);
  u10_2.resize  (nres_2); u11_2.resize (nres_2); u12_2.resize (nres_2);
  u13_2.resize  (nres_2); u14_2.resize (nres_2); u15_2.resize (nres_2);
  u16_2.resize  (nres_2); u17_2.resize (nres_2);
    
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &u1_2(j), &u2_2(j), &u3_2(j), &u4_2(j), &u5_2(j), &u6_2(j), &u7_2(j), &u8_2(j), &u9_2(j),
		  &u10_2(j), &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j), &u16_2(j), &u17_2(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_2 (3)\n");
	  exit (1);
	}
    }

  Freal_2.resize (nres_2, nres_2); Fimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_2(j, k), &Fimag_2(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_2 (4)\n");
	    exit (1);
	  }
      }

  Ereal_2.resize (nres_2, nres_2);  Eimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_2(j, k), &Eimag_2(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error reading fFile_2 (5)\n");
	    exit (1);
	  }
      }
 
  fclose (file);
  
   // ......................
  // Interpolate fFile data
  // ......................
  if (NPSI_1 != NPSI_2)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error - NPSI mismatch\n");
    }
  else
    NPSI_0 = NPSI_1;

  if (NTOR_1 != NTOR_2)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuadratic: Error - NTOR mismatch\n");
    }
  else
    NTOR_0 = NTOR_1;

  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
 
  v1_0.resize     (NPSI_0); v2_0.resize  (NPSI_0); v3_0.resize  (NPSI_0);

  mres_0.resize   (nres_0);
  u1_0.resize     (nres_0); u2_0.resize  (nres_0); u3_0.resize  (nres_0);
  u4_0.resize     (nres_0); u5_0.resize  (nres_0); u6_0.resize  (nres_0);
  u7_0.resize     (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0);
  u10_0.resize    (nres_0); u11_0.resize (nres_0); u12_0.resize (nres_0);
  u13_0.resize    (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize    (nres_0); u17_0.resize (nres_0);

  Freal_0.resize  (nres_0, nres_0); Fimag_0.resize  (nres_0, nres_0);
  Ereal_0.resize  (nres_0, nres_0); Eimag_0.resize  (nres_0, nres_0);

  if (nres_0 == nres_1)
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_1(i);
    }
  else 
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_2(i);
    }
 
  double weight1 = (time - time2) /(time1 - time2);
  double weight2 = (time - time1) /(time2 - time1);
 
  r1_0  = weight1 * r1_1  + weight2 * r1_2;
  r2_0  = weight1 * r2_1  + weight2 * r2_2;
  r3_0  = weight1 * r3_1  + weight2 * r3_2;
  r4_0  = weight1 * r4_1  + weight2 * r4_2;
  r5_0  = weight1 * r5_1  + weight2 * r5_2;
  r6_0  = weight1 * r6_1  + weight2 * r6_2;
  r7_0  = weight1 * r7_1  + weight2 * r7_2;
  r8_0  = weight1 * r8_1  + weight2 * r8_2;
  r9_0  = weight1 * r9_1  + weight2 * r9_2;
  r10_0 = weight1 * r10_1 + weight2 * r10_2;
  r11_0 = weight1 * r11_1 + weight2 * r11_2;
  r12_0 = weight1 * r12_1 + weight2 * r12_2;
  r13_0 = weight1 * r13_1 + weight2 * r13_2;

  for (int j = 0; j < NPSI_0; j++)
    {
      v1_0(j) = weight1 * v1_1(j) + weight2 * v1_2(j);
      v2_0(j) = weight1 * v2_1(j) + weight2 * v2_2(j);
      v3_0(j) = weight1 * v3_1(j) + weight2 * v3_2(j);
    }

  for (int j = 0; j < nres_0; j++)
    {
      u1_0  (j) = weight1 * u1_1  (j) + weight2 * u1_2  (j);
      u2_0  (j) = weight1 * u2_1  (j) + weight2 * u2_2  (j);
      u3_0  (j) = weight1 * u3_1  (j) + weight2 * u3_2  (j);
      u4_0  (j) = weight1 * u4_1  (j) + weight2 * u4_2  (j);
      u5_0  (j) = weight1 * u5_1  (j) + weight2 * u5_2  (j);
      u6_0  (j) = weight1 * u6_1  (j) + weight2 * u6_2  (j);
      u7_0  (j) = weight1 * u7_1  (j) + weight2 * u7_2  (j);
      u8_0  (j) = weight1 * u8_1  (j) + weight2 * u8_2  (j);
      u9_0  (j) = weight1 * u9_1  (j) + weight2 * u9_2  (j);
      u10_0 (j) = weight1 * u10_1 (j) + weight2 * u10_2 (j);
      u11_0 (j) = weight1 * u11_1 (j) + weight2 * u11_2 (j);
      u12_0 (j) = weight1 * u12_1 (j) + weight2 * u12_2 (j);
      u13_0 (j) = weight1 * u13_1 (j) + weight2 * u13_2 (j);
      u14_0 (j) = weight1 * u14_1 (j) + weight2 * u14_2 (j);
      u15_0 (j) = weight1 * u15_1 (j) + weight2 * u15_2 (j);
      u16_0 (j) = weight1 * u16_1 (j) + weight2 * u16_2 (j);
      u17_0 (j) = weight1 * u17_1 (j) + weight2 * u17_2 (j);
    }

  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      {
	Freal_0(j, k) = weight1 * Freal_1(j, k) + weight2 * Freal_2(j, k);
	Fimag_0(j, k) = weight1 * Fimag_1(j, k) + weight2 * Fimag_2(j, k);
	Ereal_0(j, k) = weight1 * Ereal_1(j, k) + weight2 * Ereal_2(j, k);
	Eimag_0(j, k) = weight1 * Eimag_1(j, k) + weight2 * Eimag_2(j, k);
      }
  // ........................
  // Write interpolated fFile
  // ........................
  file = OpenFilew (fFile);
  
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e %16.9e %16.9e %16.9e\n",
	   r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, NPSI_0, NTOR_0, nres_0, r10_0, r11_0, r12_0, r13_0);
  for (int j = 0; j < NPSI_0; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     v1_0(j), v2_0(j), v3_0(j));
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j), u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Freal_0(j, k), Fimag_0(j, k));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Ereal_0(j, k), Eimag_0(j, k));
 
  fclose (file);
  
  printf ("fFile Interpolation:\n");
  printf ("%s %11.4e\n", fFile1, weight1);
  printf ("%s %11.4e\n", fFile2, weight2);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "fFile Interpolation:\n");
  fprintf (monitor, "%s %11.4e\n", fFile1, weight1);
  fprintf (monitor, "%s %11.4e\n", fFile2, weight2);
  fclose  (monitor);
}

void Neoclassical::fFileInterpolateCubic (char* fFile1, double time1, char* fFile2, double time2, char* fFile3, double time3, char* fFile, double time)
{
  int ini; double inr;

  double          r1_1, r2_1, r3_1, r4_1, r5_1, r6_1, r7_1, r8_1, r9_1, r10_1, r11_1, r12_1, r13_1;
  int             NPSI_1, NTOR_1, nres_1;
  Array<double,1> v1_1, v2_1, v3_1;
  Array<int, 1>   mres_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1, u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1;
  Array<double,2> Freal_1, Fimag_1, Ereal_1, Eimag_1;
  
  double          r1_2, r2_2, r3_2, r4_2, r5_2, r6_2, r7_2, r8_2, r9_2, r10_2, r11_2, r12_2, r13_2;
  int             NPSI_2, NTOR_2, nres_2;
  Array<double,1> v1_2, v2_2, v3_2;
  Array<int, 1>   mres_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2, u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2;
  Array<double,2> Freal_2, Fimag_2, Ereal_2, Eimag_2;
 
  double          r1_3, r2_3, r3_3, r4_3, r5_3, r6_3, r7_3, r8_3, r9_3, r10_3, r11_3, r12_3, r13_3;
  int             NPSI_3, NTOR_3, nres_3;
  Array<double,1> v1_3, v2_3, v3_3;
  Array<int, 1>   mres_3;
  Array<double,1> u1_3, u2_3, u3_3, u4_3, u5_3, u6_3, u7_3, u8_3, u9_3, u10_3, u11_3, u12_3, u13_3, u14_3, u15_3, u16_3, u17_3;
  Array<double,2> Freal_3, Fimag_3, Ereal_3, Eimag_3;

  double          r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, r10_0, r11_0, r12_0, r13_0;
  int             NPSI_0, NTOR_0, nres_0;
  Array<double,1> v1_0, v2_0, v3_0;
  Array<int, 1>   mres_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0, u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0;
  Array<double,2> Freal_0, Fimag_0, Ereal_0, Eimag_0;

  // ................
  // Read first fFile
  // ................
  FILE* file = OpenFiler (fFile1);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_1, &r2_1, &r3_1, &r4_1, &r5_1, &r6_1, &r7_1, &r8_1, &r9_1, &NPSI_1, &NTOR_1, &nres_1, &r10_1, &r11_1, &r12_1, &r13_1) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_1 (1)\n");
      exit (1);
    }

  v1_1.resize (NPSI_1); v2_1.resize (NPSI_1); v3_1.resize (NPSI_1);

  for (int j = 0; j < NPSI_1; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_1(j), &v2_1(j), &v3_1(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_1 (2)\n");
	  exit (1);
	}
    }

  mres_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize  (nres_1); u3_1.resize  (nres_1);
  u4_1.resize   (nres_1); u5_1.resize  (nres_1); u6_1.resize  (nres_1);
  u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1);
  u10_1.resize  (nres_1); u11_1.resize (nres_1); u12_1.resize (nres_1);
  u13_1.resize  (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize (nres_1);
  
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &u1_1(j), &u2_1(j), &u3_1(j), &u4_1(j), &u5_1(j), &u6_1(j), &u7_1(j), &u8_1(j), &u9_1(j),
		  &u10_1(j), &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j), &u16_1(j), &u17_1(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_1 (3)\n");
	  exit (1);
	}
    }

  Freal_1.resize (nres_1, nres_1); Fimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_1(j, k), &Fimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_1 (4)\n");
	    exit (1);
	  }
      }
  
  Ereal_1.resize (nres_1, nres_1); Eimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_1(j, k), &Eimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_1 (5)\n");
	    exit (1);
	  }
      }
  
  fclose (file);

  // .................
  // Read second fFile
  // .................
  file = OpenFiler (fFile2);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_2, &r2_2, &r3_2, &r4_2, &r5_2, &r6_2, &r7_2, &r8_2, &r9_2, &NPSI_2, &NTOR_2, &nres_2, &r10_2, &r11_2, &r12_2, &r13_2) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_2 (1)\n");
      exit (1);
    }

  v1_2.resize (NPSI_2); v2_2.resize (NPSI_2); v3_2.resize (NPSI_2);

  for (int j = 0; j < NPSI_2; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_2(j), &v2_2(j), &v3_2(j)) != 3)
	{
	  printf ("'NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_2 (2)\n");
	  exit (1);
	}
    }

  mres_2.resize (nres_2);
  u1_2.resize   (nres_2); u2_2.resize  (nres_2); u3_2.resize  (nres_2);
  u4_2.resize   (nres_2); u5_2.resize  (nres_2); u6_2.resize  (nres_2);
  u7_2.resize   (nres_2); u8_2.resize  (nres_2); u9_2.resize  (nres_2);
  u10_2.resize  (nres_2); u11_2.resize (nres_2); u12_2.resize (nres_2);
  u13_2.resize  (nres_2); u14_2.resize (nres_2); u15_2.resize (nres_2);
  u16_2.resize  (nres_2); u17_2.resize (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &u1_2(j), &u2_2(j), &u3_2(j), &u4_2(j), &u5_2(j), &u6_2(j), &u7_2(j), &u8_2(j), &u9_2(j),
		  &u10_2(j), &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j), &u16_2(j), &u17_2(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_2 (3)\n");
	  exit (1);
	}
    }

  Freal_2.resize (nres_2, nres_2);
  Fimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_2(j, k), &Fimag_2(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_2 (4)\n");
	    exit (1);
	  }
      }

  Ereal_2.resize (nres_2, nres_2); Eimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_2(j, k), &Eimag_2(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_2 (5)\n");
	    exit (1);
	  }
      }

  fclose (file);
  
  // ................
  // Read third fFile
  // ................
  file = OpenFiler (fFile3);
  
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_3, &r2_3, &r3_3, &r4_3, &r5_3, &r6_3, &r7_3, &r8_3, &r9_3, &NPSI_3, &NTOR_3, &nres_3, &r10_3, &r11_3, &r12_3, &r13_3) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_3 (1)\n");
      exit (1);
    }

  v1_3.resize (NPSI_3); v2_3.resize (NPSI_3); v3_3.resize (NPSI_3);

  for (int j = 0; j < NPSI_3; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_3(j), &v2_3(j), &v3_3(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_3 (2)\n");
	  exit (1);
	}
    }

  mres_3.resize (nres_3);
  u1_3.resize   (nres_3); u2_3.resize  (nres_3); u3_3.resize  (nres_3);
  u4_3.resize   (nres_3); u5_3.resize  (nres_3); u6_3.resize  (nres_3);
  u7_3.resize   (nres_3); u8_3.resize  (nres_3); u9_3.resize  (nres_3);
  u10_3.resize  (nres_3); u11_3.resize (nres_3); u12_3.resize (nres_3);
  u13_3.resize  (nres_3); u14_3.resize (nres_3); u15_3.resize (nres_3);
  u16_3.resize  (nres_3); u17_3.resize (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &u1_3(j), &u2_3(j), &u3_3(j), &u4_3(j), &u5_3(j), &u6_3(j), &u7_3(j), &u8_3(j), &u9_3(j),
		  &u10_3(j), &u11_3(j), &u12_3(j), &u13_3(j), &u14_3(j), &u15_3(j), &u16_3(j), &u17_3(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_3 (3)\n");
	  exit (1);
	}
    }

  Freal_3.resize (nres_3, nres_3); Fimag_3.resize (nres_3, nres_3);
  
  for (int j = 0; j < nres_3; j++)
    for (int k = 0; k < nres_3; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_3(j, k), &Fimag_3(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_3 (4)\n");
	    exit (1);
	  }
      }

  Ereal_3.resize (nres_3, nres_3); Eimag_3.resize (nres_3, nres_3);
  
  for (int j = 0; j < nres_3; j++)
    for (int k = 0; k < nres_3; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_3(j, k), &Eimag_3(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateCubic: Error reading fFile_3 (5)\n");
	    exit (1);
	  }
      }
   
  fclose (file);
  
  // ......................
  // Interpolate fFile data
  // ......................
  if (NPSI_1 != NPSI_2 || NPSI_2 != NPSI_3)
    {
      printf ("NEOCLASSICAL::fFileInterpolateCubic: Error - NPSI mismatch\n");
    }
  else
    NPSI_0 = NPSI_1;

  if (NTOR_1 != NTOR_2 || NTOR_2 != NTOR_3)
    {
      printf ("NEOCLASSICAL::fFileInterpolateCubic: Error - NTOR mismatch\n");
    }
  else
    NTOR_0 = NTOR_1;

  if (nres_1 == nres_2 && nres_2 != nres_3)
    {
      fFileInterpolateQuadratic (fFile1, time1, fFile2, time2, fFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2)
    {
      fFileInterpolateQuadratic (fFile2, time2, fFile3, time3, fFile, time);
      return;
    }
    
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
  if (nres_3 < nres_0)
    nres_0 = nres_3;
  
  v1_0.resize     (NPSI_0); v2_0.resize  (NPSI_0); v3_0.resize  (NPSI_0);

  mres_0.resize   (nres_0);

  u1_0.resize     (nres_0); u2_0.resize  (nres_0); u3_0.resize  (nres_0);
  u4_0.resize     (nres_0); u5_0.resize  (nres_0); u6_0.resize  (nres_0);
  u7_0.resize     (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0);
  u10_0.resize    (nres_0); u11_0.resize (nres_0); u12_0.resize (nres_0);
  u13_0.resize    (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize    (nres_0); u17_0.resize (nres_0);

  Freal_0.resize  (nres_0, nres_0); Fimag_0.resize  (nres_0, nres_0);
  Ereal_0.resize  (nres_0, nres_0); Eimag_0.resize  (nres_0, nres_0);

  if (nres_0 == nres_1)
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_1(i);
    }
  else if (nres_0 == nres_2)
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_2(i);
    }
  else
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_3(i);
    }
   
  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);

  r1_0  = weight1 * r1_1  + weight2 * r1_2  + weight3 * r1_3;
  r2_0  = weight1 * r2_1  + weight2 * r2_2  + weight3 * r2_3;
  r3_0  = weight1 * r3_1  + weight2 * r3_2  + weight3 * r3_3;
  r4_0  = weight1 * r4_1  + weight2 * r4_2  + weight3 * r4_3;
  r5_0  = weight1 * r5_1  + weight2 * r5_2  + weight3 * r5_3;
  r6_0  = weight1 * r6_1  + weight2 * r6_2  + weight3 * r6_3;
  r7_0  = weight1 * r7_1  + weight2 * r7_2  + weight3 * r7_3;
  r8_0  = weight1 * r8_1  + weight2 * r8_2  + weight3 * r8_3;
  r9_0  = weight1 * r9_1  + weight2 * r9_2  + weight3 * r9_3;
  r10_0 = weight1 * r10_1 + weight2 * r10_2 + weight3 * r10_3;
  r11_0 = weight1 * r11_1 + weight2 * r11_2 + weight3 * r11_3;
  r12_0 = weight1 * r12_1 + weight2 * r12_2 + weight3 * r12_3;
  r13_0 = weight1 * r13_1 + weight2 * r13_2 + weight3 * r13_3;
  
  for (int j = 0; j < NPSI_0; j++)
    {
      v1_0(j) = weight1 * v1_1(j) + weight2 * v1_2(j) + weight3 * v1_3(j);
      v2_0(j) = weight1 * v2_1(j) + weight2 * v2_2(j) + weight3 * v2_3(j);
      v3_0(j) = weight1 * v3_1(j) + weight2 * v3_2(j) + weight3 * v3_3(j);
    }

  for (int j = 0; j < nres_0; j++)
    {
      u1_0  (j) = weight1 * u1_1  (j) + weight2 * u1_2  (j) + weight3 * u1_3  (j);
      u2_0  (j) = weight1 * u2_1  (j) + weight2 * u2_2  (j) + weight3 * u2_3  (j);
      u3_0  (j) = weight1 * u3_1  (j) + weight2 * u3_2  (j) + weight3 * u3_3  (j);
      u4_0  (j) = weight1 * u4_1  (j) + weight2 * u4_2  (j) + weight3 * u4_3  (j);
      u5_0  (j) = weight1 * u5_1  (j) + weight2 * u5_2  (j) + weight3 * u5_3  (j);
      u6_0  (j) = weight1 * u6_1  (j) + weight2 * u6_2  (j) + weight3 * u6_3  (j);
      u7_0  (j) = weight1 * u7_1  (j) + weight2 * u7_2  (j) + weight3 * u7_3  (j);
      u8_0  (j) = weight1 * u8_1  (j) + weight2 * u8_2  (j) + weight3 * u8_3  (j);
      u9_0  (j) = weight1 * u9_1  (j) + weight2 * u9_2  (j) + weight3 * u9_3  (j);
      u10_0 (j) = weight1 * u10_1 (j) + weight2 * u10_2 (j) + weight3 * u10_3 (j);
      u11_0 (j) = weight1 * u11_1 (j) + weight2 * u11_2 (j) + weight3 * u11_3 (j);
      u12_0 (j) = weight1 * u12_1 (j) + weight2 * u12_2 (j) + weight3 * u12_3 (j);
      u13_0 (j) = weight1 * u13_1 (j) + weight2 * u13_2 (j) + weight3 * u13_3 (j);
      u14_0 (j) = weight1 * u14_1 (j) + weight2 * u14_2 (j) + weight3 * u14_3 (j);
      u15_0 (j) = weight1 * u15_1 (j) + weight2 * u15_2 (j) + weight3 * u15_3 (j);
      u16_0 (j) = weight1 * u16_1 (j) + weight2 * u16_2 (j) + weight3 * u16_3 (j);
      u17_0 (j) = weight1 * u17_1 (j) + weight2 * u17_2 (j) + weight3 * u17_3 (j);
    }

  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      {
	Freal_0(j, k) = weight1 * Freal_1(j, k) + weight2 * Freal_2(j, k) + weight3 * Freal_3(j, k);
	Fimag_0(j, k) = weight1 * Fimag_1(j, k) + weight2 * Fimag_2(j, k) + weight3 * Fimag_3(j, k);
	Ereal_0(j, k) = weight1 * Ereal_1(j, k) + weight2 * Ereal_2(j, k) + weight3 * Ereal_3(j, k);
	Eimag_0(j, k) = weight1 * Eimag_1(j, k) + weight2 * Eimag_2(j, k) + weight3 * Eimag_3(j, k);
      }

  // ........................
  // Write interpolated fFile
  // ........................
  file = OpenFilew (fFile);
  
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e %16.9e %16.9e %16.9e\n",
	   r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, NPSI_0, NTOR_0, nres_0, r10_0, r11_0, r12_0, r13_0);
  for (int j = 0; j < NPSI_0; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     v1_0(j), v2_0(j), v3_0(j));
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j), u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Freal_0(j, k), Fimag_0(j, k));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Ereal_0(j, k), Eimag_0(j, k));
 
  fclose (file);
  
  printf ("fFile Interpolation:\n");
  printf ("%s %11.4e\n", fFile1, weight1);
  printf ("%s %11.4e\n", fFile2, weight2);
  printf ("%s %11.4e\n", fFile3, weight3);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "fFile Interpolation:\n");
  fprintf (monitor, "%s %11.4e\n", fFile1, weight1);
  fprintf (monitor, "%s %11.4e\n", fFile2, weight2);
  fprintf (monitor, "%s %11.4e\n", fFile3, weight3);
  fclose (monitor);
}

void Neoclassical::fFileInterpolateQuartic (char* fFile1, double time1, char* fFile2, double time2, char* fFile3, double time3,
					    char* fFile4, double time4, char* fFile, double time)
{
  int ini; double inr;

  double          r1_1, r2_1, r3_1, r4_1, r5_1, r6_1, r7_1, r8_1, r9_1, r10_1, r11_1, r12_1, r13_1;
  int             NPSI_1, NTOR_1, nres_1;
  Array<double,1> v1_1, v2_1, v3_1;
  Array<int, 1>   mres_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1, u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1;
  Array<double,2> Freal_1, Fimag_1, Ereal_1, Eimag_1;
  
  double          r1_2, r2_2, r3_2, r4_2, r5_2, r6_2, r7_2, r8_2, r9_2, r10_2, r11_2, r12_2, r13_2;
  int             NPSI_2, NTOR_2, nres_2;
  Array<double,1> v1_2, v2_2, v3_2;
  Array<int, 1>   mres_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2, u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2;
  Array<double,2> Freal_2, Fimag_2, Ereal_2, Eimag_2;
 
  double          r1_3, r2_3, r3_3, r4_3, r5_3, r6_3, r7_3, r8_3, r9_3, r10_3, r11_3, r12_3, r13_3;
  int             NPSI_3, NTOR_3, nres_3;
  Array<double,1> v1_3, v2_3, v3_3;
  Array<int, 1>   mres_3;
  Array<double,1> u1_3, u2_3, u3_3, u4_3, u5_3, u6_3, u7_3, u8_3, u9_3, u10_3, u11_3, u12_3, u13_3, u14_3, u15_3, u16_3, u17_3;
  Array<double,2> Freal_3, Fimag_3, Ereal_3, Eimag_3;
 
  double          r1_4, r2_4, r3_4, r4_4, r5_4, r6_4, r7_4, r8_4, r9_4, r10_4, r11_4, r12_4, r13_4;
  int             NPSI_4, NTOR_4, nres_4;
  Array<double,1> v1_4, v2_4, v3_4;
  Array<int, 1>   mres_4;
  Array<double,1> u1_4, u2_4, u3_4, u4_4, u5_4, u6_4, u7_4, u8_4, u9_4, u10_4, u11_4, u12_4, u13_4, u14_4, u15_4, u16_4, u17_4;
  Array<double,2> Freal_4, Fimag_4, Ereal_4, Eimag_4;

  double          r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, r10_0, r11_0, r12_0, r13_0;
  int             NPSI_0, NTOR_0, nres_0;
  Array<double,1> v1_0, v2_0, v3_0;
  Array<int, 1>   mres_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0, u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0;
  Array<double,2> Freal_0, Fimag_0, Ereal_0, Eimag_0;

  // ................
  // Read first fFile
  // ................
  FILE* file = OpenFiler (fFile1);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_1, &r2_1, &r3_1, &r4_1, &r5_1, &r6_1, &r7_1, &r8_1, &r9_1, &NPSI_1, &NTOR_1, &nres_1, &r10_1, &r11_1, &r12_1, &r13_1) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_1 (1)\n");
      exit (1);
    }

  v1_1.resize (NPSI_1); v2_1.resize (NPSI_1); v3_1.resize (NPSI_1);

  for (int j = 0; j < NPSI_1; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_1(j), &v2_1(j), &v3_1(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_1 (2)\n");
	  exit (1);
	}
    }

  mres_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize  (nres_1); u3_1.resize  (nres_1);
  u4_1.resize   (nres_1); u5_1.resize  (nres_1); u6_1.resize  (nres_1);
  u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1);
  u10_1.resize  (nres_1); u11_1.resize (nres_1); u12_1.resize (nres_1);
  u13_1.resize  (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize (nres_1);
  
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &u1_1(j), &u2_1(j), &u3_1(j), &u4_1(j), &u5_1(j), &u6_1(j), &u7_1(j), &u8_1(j), &u9_1(j),
		  &u10_1(j), &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j), &u16_1(j), &u17_1(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_1 (3)\n");
	  exit (1);
	}
    }

  Freal_1.resize (nres_1, nres_1); Fimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_1(j, k), &Fimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_1 (4)\n");
	    exit (1);
	  }
      }
  
  Ereal_1.resize (nres_1, nres_1); Eimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_1(j, k), &Eimag_1(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_1 (5)\n");
	    exit (1);
	  }
      }
  
  fclose (file);

  // .................
  // Read second fFile
  // .................
  file = OpenFiler (fFile2);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_2, &r2_2, &r3_2, &r4_2, &r5_2, &r6_2, &r7_2, &r8_2, &r9_2, &NPSI_2, &NTOR_2, &nres_2, &r10_2, &r11_2, &r12_2, &r13_2) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_2 (1)\n");
      exit (1);
    }

  v1_2.resize (NPSI_2); v2_2.resize (NPSI_2); v3_2.resize (NPSI_2);

  for (int j = 0; j < NPSI_2; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_2(j), &v2_2(j), &v3_2(j)) != 3)
	{
	  printf ("'NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_2 (2)\n");
	  exit (1);
	}
    }

  mres_2.resize (nres_2);
  u1_2.resize   (nres_2); u2_2.resize  (nres_2); u3_2.resize  (nres_2);
  u4_2.resize   (nres_2); u5_2.resize  (nres_2); u6_2.resize  (nres_2);
  u7_2.resize   (nres_2); u8_2.resize  (nres_2); u9_2.resize  (nres_2);
  u10_2.resize  (nres_2); u11_2.resize (nres_2); u12_2.resize (nres_2);
  u13_2.resize  (nres_2); u14_2.resize (nres_2); u15_2.resize (nres_2);
  u16_2.resize  (nres_2); u17_2.resize (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &u1_2(j), &u2_2(j), &u3_2(j), &u4_2(j), &u5_2(j), &u6_2(j), &u7_2(j), &u8_2(j), &u9_2(j),
		  &u10_2(j), &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j), &u16_2(j), &u17_2(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_2 (3)\n");
	  exit (1);
	}
    }

  Freal_2.resize (nres_2, nres_2); Fimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_2(j, k), &Fimag_2(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_2 (4)\n");
	    exit (1);
	  }
      }

  Ereal_2.resize (nres_2, nres_2); Eimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_2(j, k), &Eimag_2(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_2 (5)\n");
	    exit (1);
	  }
      }
  
  fclose (file);
  
  // ................
  // Read third fFile
  // ................
  file = OpenFiler (fFile3);
  
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_3, &r2_3, &r3_3, &r4_3, &r5_3, &r6_3, &r7_3, &r8_3, &r9_3, &NPSI_3, &NTOR_3, &nres_3, &r10_3, &r11_3, &r12_3, &r13_3) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_3 (1)\n");
      exit (1);
    }

  v1_3.resize (NPSI_3); v2_3.resize (NPSI_3); v3_3.resize (NPSI_3);

  for (int j = 0; j < NPSI_3; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_3(j), &v2_3(j), &v3_3(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_3 (2)\n");
	  exit (1);
	}
    }

  mres_3.resize (nres_3);
  u1_3.resize   (nres_3); u2_3.resize  (nres_3); u3_3.resize  (nres_3);
  u4_3.resize   (nres_3); u5_3.resize  (nres_3); u6_3.resize  (nres_3);
  u7_3.resize   (nres_3); u8_3.resize  (nres_3); u9_3.resize  (nres_3);
  u10_3.resize  (nres_3); u11_3.resize (nres_3); u12_3.resize (nres_3);
  u13_3.resize  (nres_3); u14_3.resize (nres_3); u15_3.resize (nres_3);
  u16_3.resize  (nres_3); u17_3.resize (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &u1_3(j), &u2_3(j), &u3_3(j), &u4_3(j), &u5_3(j), &u6_3(j), &u7_3(j), &u8_3(j), &u9_3(j),
		  &u10_3(j), &u11_3(j), &u12_3(j), &u13_3(j), &u14_3(j), &u15_3(j), &u16_3(j), &u17_3(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_3 (3)\n");
	  exit (1);
	}
    }

  Freal_3.resize (nres_3, nres_3); Fimag_3.resize (nres_3, nres_3);
  
  for (int j = 0; j < nres_3; j++)
    for (int k = 0; k < nres_3; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_3(j, k), &Fimag_3(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_3 (4)\n");
	    exit (1);
	  }
      }

  Ereal_3.resize (nres_3, nres_3);
  Eimag_3.resize (nres_3, nres_3);
  
  for (int j = 0; j < nres_3; j++)
    for (int k = 0; k < nres_3; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_3(j, k), &Eimag_3(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_3 (5)\n");
	    exit (1);
	  }
      }
  
  fclose (file);

  // .................
  // Read fourth fFile
  // .................
  file = OpenFiler (fFile4);
  
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf",
	      &r1_4, &r2_4, &r3_4, &r4_4, &r5_4, &r6_4, &r7_4, &r8_4, &r9_4, &NPSI_4, &NTOR_4, &nres_4, &r10_4, &r11_4, &r12_4, &r13_4) != 16)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_4 (1)\n");
      exit (1);
    }

  v1_4.resize (NPSI_4); v2_4.resize (NPSI_4); v3_4.resize (NPSI_4);

  for (int j = 0; j < NPSI_4; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_4(j), &v2_4(j), &v3_4(j)) != 3)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_4 (2)\n");
	  exit (1);
	}
    }

  mres_4.resize (nres_4);
  u1_4.resize   (nres_4); u2_4.resize  (nres_4); u3_4.resize  (nres_4);
  u4_4.resize   (nres_4); u5_4.resize  (nres_4); u6_4.resize  (nres_4);
  u7_4.resize   (nres_4); u8_4.resize  (nres_4); u9_4.resize  (nres_4);
  u10_4.resize  (nres_4); u11_4.resize (nres_4); u12_4.resize (nres_4);
  u13_4.resize  (nres_4); u14_4.resize (nres_4); u15_4.resize (nres_4);
  u16_4.resize  (nres_4); u17_4.resize (nres_4);
  
  for (int j = 0; j < nres_4; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_4(j), &u1_4(j), &u2_4(j), &u3_4(j), &u4_4(j), &u5_4(j), &u6_4(j), &u7_4(j), &u8_4(j), &u9_4(j),
		  &u10_4(j), &u11_4(j), &u12_4(j), &u13_4(j), &u14_4(j), &u15_4(j), &u16_4(j), &u17_4(j)) != 18)
	{
	  printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_4 (3)\n");
	  exit (1);
	}
    }

  Freal_4.resize (nres_4, nres_4); Fimag_4.resize (nres_4, nres_4);
  
  for (int j = 0; j < nres_4; j++)
    for (int k = 0; k < nres_4; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_4(j, k), &Fimag_4(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_4 (4)\n");
	    exit (1);
	  }
      }

  Ereal_4.resize (nres_4, nres_4); Eimag_4.resize (nres_4, nres_4);
  
  for (int j = 0; j < nres_4; j++)
    for (int k = 0; k < nres_4; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_4(j, k), &Eimag_4(j, k)) != 4)
	  {
	    printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error reading fFile_4 (5)\n");
	    exit (1);
	  }
      }
 
  fclose (file);
  
  // ......................
  // Interpolate fFile data
  // ......................
  if (NPSI_1 != NPSI_2 || NPSI_2 != NPSI_3 || NPSI_3 != NPSI_4)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error - NPSI mismatch\n");
    }
  else
    NPSI_0 = NPSI_1;

  if (NTOR_1 != NTOR_2 || NTOR_2 != NTOR_3 || NTOR_3 != NTOR_4)
    {
      printf ("NEOCLASSICAL::fFileInterpolateQuartic: Error - NTOR mismatch\n");
    }
  else
    NTOR_0 = NTOR_1;

  if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 != nres_2)
    {
      fFileInterpolateQuadratic (fFile2, time2, fFile3, time3, fFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 == nres_2 && nres_4 != nres_2)
    {
      fFileInterpolateCubic (fFile1, time1, fFile2, time2, fFile3, time3, fFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 == nres_2)
    {
      fFileInterpolateCubic (fFile2, time2, fFile3, time3, fFile4, time4, fFile, time);
      return;
    }
  
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
  if (nres_3 < nres_0)
    nres_0 = nres_3;
  if (nres_4 < nres_0)
    nres_0 = nres_4;
  
  v1_0.resize     (NPSI_0); v2_0.resize  (NPSI_0); v3_0.resize  (NPSI_0);

  mres_0.resize   (nres_0);
  u1_0.resize     (nres_0); u2_0.resize  (nres_0); u3_0.resize  (nres_0);
  u4_0.resize     (nres_0); u5_0.resize  (nres_0); u6_0.resize  (nres_0);
  u7_0.resize     (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0);
  u10_0.resize    (nres_0); u11_0.resize (nres_0); u12_0.resize (nres_0);
  u13_0.resize    (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize    (nres_0); u17_0.resize (nres_0);

  Freal_0.resize  (nres_0, nres_0); Fimag_0.resize  (nres_0, nres_0);
  Ereal_0.resize  (nres_0, nres_0); Eimag_0.resize  (nres_0, nres_0);

  if (nres_0 == nres_1)
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_1(i);
    }
  else if (nres_0 == nres_2)
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_2(i);
    }
  else if (nres_0 == nres_3)
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_3(i);
    }
  else
    {
      for (int i = 0; i < nres_0; i++)
	mres_0(i) = mres_4(i);
    }

  double weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4);
  double weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4);
  double weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4);
  double weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3);

  r1_0  = weight1 * r1_1  + weight2 * r1_2  + weight3 * r1_3  + weight4 * r1_4;
  r2_0  = weight1 * r2_1  + weight2 * r2_2  + weight3 * r2_3  + weight4 * r2_4;
  r3_0  = weight1 * r3_1  + weight2 * r3_2  + weight3 * r3_3  + weight4 * r3_4;
  r4_0  = weight1 * r4_1  + weight2 * r4_2  + weight3 * r4_3  + weight4 * r4_4;
  r5_0  = weight1 * r5_1  + weight2 * r5_2  + weight3 * r5_3  + weight4 * r5_4;
  r6_0  = weight1 * r6_1  + weight2 * r6_2  + weight3 * r6_3  + weight4 * r6_4;
  r7_0  = weight1 * r7_1  + weight2 * r7_2  + weight3 * r7_3  + weight4 * r7_4;
  r8_0  = weight1 * r8_1  + weight2 * r8_2  + weight3 * r8_3  + weight4 * r8_4;
  r9_0  = weight1 * r9_1  + weight2 * r9_2  + weight3 * r9_3  + weight4 * r9_4;
  r10_0 = weight1 * r10_1 + weight2 * r10_2 + weight3 * r10_3 + weight4 * r10_4;
  r11_0 = weight1 * r11_1 + weight2 * r11_2 + weight3 * r11_3 + weight4 * r11_4;
  r12_0 = weight1 * r12_1 + weight2 * r12_2 + weight3 * r12_3 + weight4 * r12_4;
  r13_0 = weight1 * r13_1 + weight2 * r13_2 + weight3 * r13_3 + weight4 * r13_4;  
 
  for (int j = 0; j < NPSI_0; j++)
    {
      v1_0(j) = weight1 * v1_1(j) + weight2 * v1_2(j) + weight3 * v1_3(j) + weight4 * v1_4(j);
      v2_0(j) = weight1 * v2_1(j) + weight2 * v2_2(j) + weight3 * v2_3(j) + weight4 * v2_4(j);
      v3_0(j) = weight1 * v3_1(j) + weight2 * v3_2(j) + weight3 * v3_3(j) + weight4 * v3_4(j);
    }

  for (int j = 0; j < nres_0; j++)
    {
      u1_0  (j) = weight1 * u1_1  (j) + weight2 * u1_2  (j) + weight3 * u1_3  (j) + weight4 * u1_4  (j);
      u2_0  (j) = weight1 * u2_1  (j) + weight2 * u2_2  (j) + weight3 * u2_3  (j) + weight4 * u2_4  (j);
      u3_0  (j) = weight1 * u3_1  (j) + weight2 * u3_2  (j) + weight3 * u3_3  (j) + weight4 * u3_4  (j);
      u4_0  (j) = weight1 * u4_1  (j) + weight2 * u4_2  (j) + weight3 * u4_3  (j) + weight4 * u4_4  (j);
      u5_0  (j) = weight1 * u5_1  (j) + weight2 * u5_2  (j) + weight3 * u5_3  (j) + weight4 * u5_4  (j);
      u6_0  (j) = weight1 * u6_1  (j) + weight2 * u6_2  (j) + weight3 * u6_3  (j) + weight4 * u6_4  (j);
      u7_0  (j) = weight1 * u7_1  (j) + weight2 * u7_2  (j) + weight3 * u7_3  (j) + weight4 * u7_4  (j);
      u8_0  (j) = weight1 * u8_1  (j) + weight2 * u8_2  (j) + weight3 * u8_3  (j) + weight4 * u8_4  (j);
      u9_0  (j) = weight1 * u9_1  (j) + weight2 * u9_2  (j) + weight3 * u9_3  (j) + weight4 * u9_4  (j);
      u10_0 (j) = weight1 * u10_1 (j) + weight2 * u10_2 (j) + weight3 * u10_3 (j) + weight4 * u10_4 (j);
      u11_0 (j) = weight1 * u11_1 (j) + weight2 * u11_2 (j) + weight3 * u11_3 (j) + weight4 * u11_4 (j);
      u12_0 (j) = weight1 * u12_1 (j) + weight2 * u12_2 (j) + weight3 * u12_3 (j) + weight4 * u12_4 (j);
      u13_0 (j) = weight1 * u13_1 (j) + weight2 * u13_2 (j) + weight3 * u13_3 (j) + weight4 * u13_4 (j);
      u14_0 (j) = weight1 * u14_1 (j) + weight2 * u14_2 (j) + weight3 * u14_3 (j) + weight4 * u14_4 (j);
      u15_0 (j) = weight1 * u15_1 (j) + weight2 * u15_2 (j) + weight3 * u15_3 (j) + weight4 * u15_4 (j);
      u16_0 (j) = weight1 * u16_1 (j) + weight2 * u16_2 (j) + weight3 * u16_3 (j) + weight4 * u16_4 (j);
      u17_0 (j) = weight1 * u17_1 (j) + weight2 * u17_2 (j) + weight3 * u17_3 (j) + weight4 * u17_4 (j);
    }

  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      {
	Freal_0(j, k) = weight1 * Freal_1(j, k) + weight2 * Freal_2(j, k) + weight3 * Freal_3(j, k) + weight4 * Freal_4(j, k);
	Fimag_0(j, k) = weight1 * Fimag_1(j, k) + weight2 * Fimag_2(j, k) + weight3 * Fimag_3(j, k) + weight4 * Fimag_4(j, k);
	Ereal_0(j, k) = weight1 * Ereal_1(j, k) + weight2 * Ereal_2(j, k) + weight3 * Ereal_3(j, k) + weight4 * Ereal_4(j, k);
	Eimag_0(j, k) = weight1 * Eimag_1(j, k) + weight2 * Eimag_2(j, k) + weight3 * Eimag_3(j, k) + weight4 * Eimag_4(j, k);
      }
  
  // ........................
  // Write interpolated fFile
  // ........................
  file = OpenFilew (fFile);
  
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e %16.9e %16.9e %16.9e\n",
	   r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, r9_0, NPSI_0, NTOR_0, nres_0, r10_0, r11_0, r12_0, r13_0);
  for (int j = 0; j < NPSI_0; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     v1_0(j), v2_0(j), v3_0(j));
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j), u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Freal_0(j, k), Fimag_0(j, k));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Ereal_0(j, k), Eimag_0(j, k));
 
  fclose (file);
  
  printf ("fFile Interpolation:\n");
  printf ("%s %11.4e\n", fFile1, weight1);
  printf ("%s %11.4e\n", fFile2, weight2);
  printf ("%s %11.4e\n", fFile3, weight3);
  printf ("%s %11.4e\n", fFile4, weight4);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "fFile Interpolation:\n");
  fprintf (monitor, "%s %11.4e\n", fFile1, weight1);
  fprintf (monitor, "%s %11.4e\n", fFile2, weight2);
  fprintf (monitor, "%s %11.4e\n", fFile3, weight3);
  fprintf (monitor, "%s %11.4e\n", fFile4, weight4);
  fclose  (monitor);
}
  
