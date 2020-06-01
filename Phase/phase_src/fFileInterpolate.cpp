// fFileInterpolate.h

#include "Phase.h"

// ###############################
// Functions to interpolate fFiles
// ###############################
void Phase::fFileInterpolate (char* fFile1, double time1, char* fFile2, double time2, char* fFile3, double time3, char* fFile, double time)
{
  int ini; double inr;

  double          r1_1, r2_1, r3_1, r4_1, r5_1, r6_1, r7_1;
  int             NPSI_1, NTOR_1, nres_1;
  Array<double,1> v1_1, v2_1, v3_1;
  Array<int, 1>   mres_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1;
  Array<double,2> Freal_1, Fimag_1, Ereal_1, Eimag_1;
  Array<double,1> EIreal_1, EIimag_1, EOreal_1, EOimag_1;
  
  double          r1_2, r2_2, r3_2, r4_2, r5_2, r6_2, r7_2;
  int             NPSI_2, NTOR_2, nres_2;
  Array<double,1> v1_2, v2_2, v3_2;
  Array<int, 1>   mres_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2;
  Array<double,2> Freal_2, Fimag_2, Ereal_2, Eimag_2;
  Array<double,1> EIreal_2, EIimag_2, EOreal_2, EOimag_2;
 
  double          r1_3, r2_3, r3_3, r4_3, r5_3, r6_3, r7_3;
  int             NPSI_3, NTOR_3, nres_3;
  Array<double,1> v1_3, v2_3, v3_3;
  Array<int, 1>   mres_3;
  Array<double,1> u1_3, u2_3, u3_3, u4_3, u5_3, u6_3, u7_3, u8_3, u9_3, u10_3;
  Array<double,2> Freal_3, Fimag_3, Ereal_3, Eimag_3;
  Array<double,1> EIreal_3, EIimag_3, EOreal_3, EOimag_3;

  double          r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0;
  int             NPSI_0, NTOR_0, nres_0;
  Array<double,1> v1_0, v2_0, v3_0;
  Array<int, 1>   mres_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0;
  Array<double,2> Freal_0, Fimag_0, Ereal_0, Eimag_0;
  Array<double,1> EIreal_0, EIimag_0, EOreal_0, EOimag_0;

  // ................
  // Read first fFile
  // ................
  FILE* file = OpenFiler (fFile1);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %d %d %d",
	      &r1_1, &r2_1, &r3_1, &r4_1, &r5_1, &r6_1, &r7_1, &NPSI_1, &NTOR_1, &nres_1) != 10)
    {
      printf ("PHASE::Error reading fFile_1 (1)\n");
      exit (1);
    }

  v1_1.resize (NPSI_1);
  v2_1.resize (NPSI_1);
  v3_1.resize (NPSI_1);

  for (int j = 0; j < NPSI_1; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_1(j), &v2_1(j), &v3_1(j)) != 3)
	{
	  printf ("PHASE: Error reading fFile_1 (2)\n");
	  exit (1);
	}
    }

  mres_1.resize (nres_1);
  u1_1.resize   (nres_1);
  u2_1.resize   (nres_1);
  u3_1.resize   (nres_1);
  u4_1.resize   (nres_1);
  u5_1.resize   (nres_1);
  u6_1.resize   (nres_1);
  u7_1.resize   (nres_1);
  u8_1.resize   (nres_1);
  u9_1.resize   (nres_1);
  u10_1.resize  (nres_1);

  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &u1_1(j), &u2_1(j), &u3_1(j), &u4_1(j), &u5_1(j), &u6_1(j), &u7_1(j), &u8_1(j), &u9_1(j), &u10_1(j)) != 11)
	{
	  printf ("PHASE: Error reading fFile_1 (3)\n");
	  exit (1);
	}
    }

  Freal_1.resize (nres_1, nres_1);
  Fimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_1(j, k), &Fimag_1(j, k)) != 4)
	  {
	    printf ("PHASE::Error reading fFile_1 (4)\n");
	    exit (1);
	  }
      }
  
  Ereal_1.resize (nres_1, nres_1);
  Eimag_1.resize (nres_1, nres_1);

  for (int j = 0; j < nres_1; j++)
    for (int k = 0; k < nres_1; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_1(j, k), &Eimag_1(j, k)) != 4)
	  {
	    printf ("PHASE::Error reading fFile_1 (5)\n");
	    exit (1);
	  }
      }
  
  EIreal_1.resize (nres_1);
  EIimag_1.resize (nres_1);
  EOreal_1.resize (nres_1);
  EOimag_1.resize (nres_1);

  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf", &ini, &EIreal_1(j), &EIimag_1(j), &EOreal_1(j), &EOimag_1(j)) != 5)
	{
	  printf ("PHASE::Error reading fFile (6)\n");
	  exit (1);
	}
    }
  
  fclose (file);

  // .................
  // Read second fFile
  // .................
  file = OpenFiler (fFile2);
 
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %d %d %d",
	      &r1_2, &r2_2, &r3_2, &r4_2, &r5_2, &r6_2, &r7_2, &NPSI_2, &NTOR_2, &nres_2) != 10)
    {
      printf ("PHASE::Error reading fFile_2 (1)\n");
      exit (1);
    }

  v1_2.resize (NPSI_2);
  v2_2.resize (NPSI_2);
  v3_2.resize (NPSI_2);

  for (int j = 0; j < NPSI_2; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_2(j), &v2_2(j), &v3_2(j)) != 3)
	{
	  printf ("'PHASE: Error reading fFile_2 (2)\n");
	  exit (1);
	}
    }

  mres_2.resize (nres_2);
  u1_2.resize   (nres_2);
  u2_2.resize   (nres_2);
  u3_2.resize   (nres_2);
  u4_2.resize   (nres_2);
  u5_2.resize   (nres_2);
  u6_2.resize   (nres_2);
  u7_2.resize   (nres_2);
  u8_2.resize   (nres_2);
  u9_2.resize   (nres_2);
  u10_2.resize  (nres_2);
    
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &u1_2(j), &u2_2(j), &u3_2(j), &u4_2(j), &u5_2(j), &u6_2(j), &u7_2(j), &u8_2(j), &u9_2(j), &u10_2(j)) != 11)
	{
	  printf ("PHASE: Error reading fFile_2 (3)\n");
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
	    printf ("PHASE::Error reading fFile_2 (4)\n");
	    exit (1);
	  }
      }

  Ereal_2.resize (nres_2, nres_2);
  Eimag_2.resize (nres_2, nres_2);
  
  for (int j = 0; j < nres_2; j++)
    for (int k = 0; k < nres_2; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Ereal_2(j, k), &Eimag_2(j, k)) != 4)
	  {
	    printf ("PHASE::Error reading fFile_2 (5)\n");
	    exit (1);
	  }
      }
  
  EIreal_2.resize (nres_2);
  EIimag_2.resize (nres_2);
  EOreal_2.resize (nres_2);
  EOimag_2.resize (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf", &ini, &EIreal_2(j), &EIimag_2(j), &EOreal_2(j), &EOimag_2(j)) != 5)
	{
	  printf ("PHASE::Error reading fFile_2 (6)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ................
  // Read third fFile
  // ................
  file = OpenFiler (fFile3);
  
  if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %d %d %d",
	      &r1_3, &r2_3, &r3_3, &r4_3, &r5_3, &r6_3, &r7_3, &NPSI_3, &NTOR_3, &nres_3) != 10)
    {
      printf ("PHASE::Error reading fFile_3 (1)\n");
      exit (1);
    }

  v1_3.resize (NPSI_3);
  v2_3.resize (NPSI_3);
  v3_3.resize (NPSI_3);

  for (int j = 0; j < NPSI_3; j++)
    {
      if (fscanf (file, "%lf %lf %lf", &v1_3(j), &v2_3(j), &v3_3(j)) != 3)
	{
	  printf ("PHASE: Error reading fFile_3 (2)\n");
	  exit (1);
	}
    }

  mres_3.resize (nres_3);
  u1_3.resize   (nres_3);
  u2_3.resize   (nres_3);
  u3_3.resize   (nres_3);
  u4_3.resize   (nres_3);
  u5_3.resize   (nres_3);
  u6_3.resize   (nres_3);
  u7_3.resize   (nres_3);
  u8_3.resize   (nres_3);
  u9_3.resize   (nres_3);
  u10_3.resize  (nres_3);
    
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &u1_3(j), &u2_3(j), &u3_3(j), &u4_3(j), &u5_3(j), &u6_3(j), &u7_3(j), &u8_3(j), &u9_3(j), &u10_3(j)) != 11)
	{
	  printf ("PHASE: Error reading fFile_3 (3)\n");
	  exit (1);
	}
    }

  Freal_3.resize (nres_3, nres_3);
  Fimag_3.resize (nres_3, nres_3);
  
  for (int j = 0; j < nres_3; j++)
    for (int k = 0; k < nres_3; k++)
      {
	if (fscanf (file, "%d %d %lf %lf", &ini, &ini, &Freal_3(j, k), &Fimag_3(j, k)) != 4)
	  {
	    printf ("PHASE::Error reading fFile_3 (4)\n");
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
	    printf ("PHASE::Error reading fFile_3 (5)\n");
	    exit (1);
	  }
      }
  
  EIreal_3.resize (nres_3);
  EIimag_3.resize (nres_3);
  EOreal_3.resize (nres_3);
  EOimag_3.resize (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %lf %lf %lf %lf", &ini, &EIreal_3(j), &EIimag_3(j), &EOreal_3(j), &EOimag_3(j)) != 5)
	{
	  printf ("PHASE::Error reading fFile_3 (6)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ......................
  // Interpolate fFile data
  // ......................
  if (NPSI_1 != NPSI_2 || NPSI_1 != NPSI_3)
    {
      printf ("PHASE: Error - NPSI mismatch\n");
    }
  else
    NPSI_0 = NPSI_1;

  if (NTOR_1 != NTOR_2 || NTOR_1 != NTOR_3)
    {
      printf ("PHASE: Error - NTOR mismatch\n");
    }
  else
    NTOR_0 = NTOR_1;
  
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
  if (nres_3 < nres_0)
    nres_0 = nres_3;
  
  v1_0.resize     (NPSI_0);
  v2_0.resize     (NPSI_0);
  v3_0.resize     (NPSI_0);
  mres_0.resize   (nres_0);
  u1_0.resize     (nres_0);
  u2_0.resize     (nres_0);
  u3_0.resize     (nres_0);
  u4_0.resize     (nres_0);
  u5_0.resize     (nres_0);
  u6_0.resize     (nres_0);
  u7_0.resize     (nres_0);
  u8_0.resize     (nres_0);
  u9_0.resize     (nres_0);
  u10_0.resize    (nres_0);
  Freal_0.resize  (nres_0, nres_0);
  Fimag_0.resize  (nres_0, nres_0);
  Ereal_0.resize  (nres_0, nres_0);
  Eimag_0.resize  (nres_0, nres_0);
  EIreal_0.resize (nres_0);
  EIimag_0.resize (nres_0);
  EOreal_0.resize (nres_0);
  EOimag_0.resize (nres_0);

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

  r1_0 = weight1 * r1_1 + weight2 * r1_2 + weight3 * r1_3;
  r2_0 = weight1 * r2_1 + weight2 * r2_2 + weight3 * r2_3;
  r3_0 = weight1 * r3_1 + weight2 * r3_2 + weight3 * r3_3;
  r4_0 = weight1 * r4_1 + weight2 * r4_2 + weight3 * r4_3;
  r5_0 = weight1 * r5_1 + weight2 * r5_2 + weight3 * r5_3;
  r6_0 = weight1 * r6_1 + weight2 * r6_2 + weight3 * r6_3;
  r7_0 = weight1 * r7_1 + weight2 * r7_2 + weight3 * r7_3;

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
    }

  for (int j = 0; j < nres_0; j++)
    {
      EIreal_0(j) = weight1 * EIreal_1(j) + weight2 * EIreal_2(j) + weight3 * EIreal_3(j);
      EIimag_0(j) = weight1 * EIimag_1(j) + weight2 * EIimag_2(j) + weight3 * EIimag_3(j);
      EOreal_0(j) = weight1 * EOreal_1(j) + weight2 * EOreal_2(j) + weight3 * EOreal_3(j);
      EOimag_0(j) = weight1 * EOimag_1(j) + weight2 * EOimag_2(j) + weight3 * EOimag_3(j);
      
      for (int k = 0; k < nres_0; k++)
	{
	  Freal_0(j, k) = weight1 * Freal_1(j, k) + weight2 * Freal_2(j, k) + weight3 * Freal_3(j, k);
	  Fimag_0(j, k) = weight1 * Fimag_1(j, k) + weight2 * Fimag_2(j, k) + weight3 * Fimag_3(j, k);
	  Ereal_0(j, k) = weight1 * Ereal_1(j, k) + weight2 * Ereal_2(j, k) + weight3 * Ereal_3(j, k);
	  Eimag_0(j, k) = weight1 * Eimag_1(j, k) + weight2 * Eimag_2(j, k) + weight3 * Eimag_3(j, k);
	}
    }

  // ........................
  // Write interpolated fFile
  // ........................
  file = OpenFilew (fFile);
  
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d\n",
	   r1_0, r2_0, r3_0, r4_0, r5_0, r6_0, r7_0, NPSI_0, NTOR_0, nres_0);
  for (int j = 0; j < NPSI_0; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     v1_0(j), v2_0(j), v3_0(j));
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Freal_0(j, k), Fimag_0(j, k));
  for (int j = 0; j < nres_0; j++)
    for (int k = 0; k < nres_0; k++)
      fprintf (file, "%d %d %16.9e %16.9e\n",
	       mres_0(j), mres_0(k), Ereal_0(j, k), Eimag_0(j, k));
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), EIreal_0(j), EIimag_0(j), EOreal_0(j), EOimag_0(j));

  fclose (file);
  
  printf ("fFile Interpolation:\n");
  printf ("%s %11.4e\n", fFile1, weight1);
  printf ("%s %11.4e\n", fFile2, weight2);
  printf ("%s %11.4e\n", fFile3, weight3);
}
  
 void Phase::fFileInterp (vector<string> fFileName, vector<double> fFileTime, int fFileNumber, double time)
{
  int    index;
  double _time;
  
  if (time < fFileTime[0])
    {
      index = 0;
      _time = fFileTime[0];
    }
  else if (time >= fFileTime[fFileNumber-1])
    {
      index = fFileNumber - 3;
      _time = fFileTime[fFileNumber-1];
    }
  else
    {
      for (int i = 0; i < fFileNumber-1; i++)
	if (time >= fFileTime[i] && time < fFileTime[i+1])
	  {
	    index = i;
	    _time = time;
	    
	    if (index > fFileNumber-3)
	      index = fFileNumber - 3;
	  }
    }

  char* fFile = "../Flux/fFile";
  char* file1 = (char*) fFileName[index  ].c_str();
  char* file2 = (char*) fFileName[index+1].c_str();
  char* file3 = (char*) fFileName[index+2].c_str();

  fFileInterpolate (file1, fFileTime[index], file2, fFileTime[index+1], file3, fFileTime[index+2], fFile, _time);
}

