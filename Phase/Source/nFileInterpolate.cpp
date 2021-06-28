// nFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Phase:: nFileInterp               (vector<string> nFileName,   vector<double> nFileTime,   int nFileNumber, double time)
// void Phase:: nFileInterpolateLinear    (char* nFile1, double time1, char* nFile,  double time)
// void Phase:: nFileInterpolateQuadratic (char* nFile1, double time1, char* nFile2, double time2, char* nFile,     double time)
// void Phase:: nFileInterpolateCubic     (char* nFile1, double time1, char* nFile2, double time2, char* nFile3,    double time3, char* nFile, double time)
// void Phase:: nFileInterpolateQuartic   (char* nFile1, double time1, char* nFile2, double time2, char* nFile3,    double time3,
//				           char* nFile4, double time4, char* nFile,  double time)

#include "Phase.h"

// ###############################
// Functions to interpolate nFiles
// ###############################
void Phase::nFileInterp (vector<string> nFileName, vector<double> nFileTime, int nFileNumber, double time)
{
  if (nFileNumber < 1)
    {
      printf ("PHASE::nFileInterp: Error - nFileNumber must be greater than zero\n");
      exit (1);
    }
  else if (nFileNumber == 1)
    {
      char* nFile = "Inputs/nFile";
      char* file1 = (char*) nFileName[0].c_str();

      nFileInterpolateLinear (file1, nFileTime[0], nFile, time);
    }
  else if (nFileNumber == 2)
    {
      char* nFile = "Inputs/nFile";
      char* file1 = (char*) nFileName[0].c_str();
      char* file2 = (char*) nFileName[1].c_str();

      nFileInterpolateQuadratic (file1, nFileTime[0], file2, nFileTime[1], nFile, time);
    }
  else if (NATS)
    {
      int index;

      if (time < nFileTime[0])
	index = 0;
      else if (time >= nFileTime[nFileNumber-1])
	index = nFileNumber - 2;
      else
	{
	  for (int i = 0; i < nFileNumber - 1; i++)
	    if (time >= nFileTime[i] && time < nFileTime[i+1])
	      index = i;
	}
      
      char* nFile = "Inputs/nFile";
      char* file1 = (char*) nFileName[index  ].c_str();
      char* file2 = (char*) nFileName[index+1].c_str();

      nFileInterpolateQuadratic (file1, nFileTime[index], file2, nFileTime[index+1], nFile, time);
    }
  else if (nFileNumber == 3)
    {
      char* nFile = "Inputs/nFile";
      char* file1 = (char*) nFileName[0].c_str();
      char* file2 = (char*) nFileName[1].c_str();
      char* file3 = (char*) nFileName[2].c_str();

      nFileInterpolateCubic (file1, nFileTime[0], file2, nFileTime[1], file3, nFileTime[2], nFile, time);
    }
  else if (nFileNumber == 4)
    {
      char* nFile = "Inputs/nFile";
      char* file1 = (char*) nFileName[0].c_str();
      char* file2 = (char*) nFileName[1].c_str();
      char* file3 = (char*) nFileName[2].c_str();
      char* file4 = (char*) nFileName[3].c_str();

      nFileInterpolateQuartic (file1, nFileTime[0], file2, nFileTime[1], file3, nFileTime[2],
			       file4, nFileTime[3], nFile, time);
    }
  else
    {
      int index, cntrl;

      if (time < nFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= nFileTime[nFileNumber-1])
	{
	  index = nFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < nFileNumber-1; i++)
	    if (time >= nFileTime[i] && time < nFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == nFileNumber-2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* nFile = "Inputs/nFile";
	  char* file1 = (char*) nFileName[index-1].c_str();
	  char* file2 = (char*) nFileName[index  ].c_str();
	  char* file3 = (char*) nFileName[index+1].c_str();
	  char* file4 = (char*) nFileName[index+2].c_str();
	  
	  nFileInterpolateQuartic (file1, nFileTime[index-1], file2, nFileTime[index], file3, nFileTime[index+1],
				   file4, nFileTime[index+2], nFile, time);
	}
      else if (cntrl == 2)
	{
	  char* nFile = "Inputs/nFile";
	  char* file1 = (char*) nFileName[index  ].c_str();
	  char* file2 = (char*) nFileName[index+1].c_str();
	  char* file3 = (char*) nFileName[index+2].c_str();
	  
	  nFileInterpolateCubic (file1, nFileTime[index], file2, nFileTime[index+1], file3, nFileTime[index+2], nFile, time);
	}
      else if (cntrl == 3)
	{
	  char* nFile = "Inputs/nFile";
	  char* file1 = (char*) nFileName[index-1].c_str();
	  char* file2 = (char*) nFileName[index  ].c_str();
	  char* file3 = (char*) nFileName[index+1].c_str();
	  
	  nFileInterpolateCubic (file1, nFileTime[index-1], file2, nFileTime[index], file3, nFileTime[index+1], nFile, time);
	}
    }
}

void Phase::nFileInterpolateLinear (char* nFile1, double time1, char* nFile, double time)
{
  int             nres_1;
  double          taua_1, P0_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1,  u2_1,  u3_1,  u4_1,  u5_1,  u6_1,  u7_1,  u8_1,  u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;
  Array<double,1> u20_1, u21_1, u22_1, u23_1, u24_1, u25_1, u26_1, u27_1, u28_1;
  Array<double,1> u29_1, u30_1, u31_1, u32_1, u33_1, u34_1, u35_1, u36_1, u37_1;
  Array<double,1> u38_1, u39_1, u40_1;
       
  int             nres_0;
  double          taua_0, P0_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0,  u2_0,  u3_0,  u4_0,  u5_0,  u6_0,  u7_0,  u8_0,  u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
  Array<double,1> u20_0, u21_0, u22_0, u23_0, u24_0, u25_0, u26_0, u27_0, u28_0;
  Array<double,1> u29_0, u30_0, u31_0, u32_0, u33_0, u34_0, u35_0, u36_0, u37_0;
  Array<double,1> u38_0, u39_0, u40_0;
 
  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf %lf", &nres_1, &taua_1, &P0_1) != 3)
    {
      printf ("PHASE::nFileInterpolateLinear: Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1); ntor_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize   (nres_1); u3_1.resize  (nres_1); u4_1.resize  (nres_1); u5_1.resize  (nres_1);
  u6_1.resize   (nres_1); u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1); u10_1.resize (nres_1);
  u11_1.resize  (nres_1); u12_1.resize  (nres_1); u13_1.resize (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize  (nres_1); u18_1.resize (nres_1); u19_1.resize (nres_1); u20_1.resize (nres_1);
  u21_1.resize  (nres_1); u22_1.resize  (nres_1); u23_1.resize (nres_1); u24_1.resize (nres_1); u25_1.resize (nres_1);
  u26_1.resize  (nres_1); u27_1.resize  (nres_1); u28_1.resize (nres_1); u29_1.resize (nres_1); u30_1.resize (nres_1);
  u31_1.resize  (nres_1); u32_1.resize  (nres_1); u33_1.resize (nres_1); u34_1.resize (nres_1); u35_1.resize (nres_1);
  u36_1.resize  (nres_1); u37_1.resize  (nres_1); u38_1.resize (nres_1); u39_1.resize (nres_1); u40_1.resize (nres_1);
  
   for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j), &u20_1(j),
		  &u21_1(j), &u22_1(j), &u23_1(j), &u25_1(j), &u25_1(j),
		  &u26_1(j), &u27_1(j), &u28_1(j), &u29_1(j), &u30_1(j),
		  &u31_1(j), &u32_1(j), &u33_1(j), &u34_1(j), &u35_1(j),
		  &u36_1(j), &u37_1(j), &u38_1(j), &u39_1(j), &u40_1(j)) != 42)
	{
	  printf ("PHASE::nFileInterpolateLinear: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);

  // ......................
  // Interpolate nFile data
  // ......................
  nres_0 = nres_1;

  mres_0.resize (nres_0); ntor_0.resize (nres_0);
  u1_0.resize   (nres_0); u2_0.resize   (nres_0); u3_0.resize  (nres_0); u9_0.resize  (nres_0); u10_0.resize (nres_0);
  u11_0.resize  (nres_0); u12_0.resize  (nres_0); u13_0.resize (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize  (nres_0); u17_0.resize  (nres_0); u18_0.resize (nres_0); u19_0.resize (nres_0); u20_0.resize (nres_0);
  u21_0.resize  (nres_0); u22_0.resize  (nres_0); u23_0.resize (nres_0); u24_0.resize (nres_0); u25_0.resize (nres_0);
  u26_0.resize  (nres_0); u27_0.resize  (nres_0); u28_0.resize (nres_0); u29_0.resize (nres_0); u30_0.resize (nres_0);
  u31_0.resize  (nres_0); u32_0.resize  (nres_0); u33_0.resize (nres_0); u34_0.resize (nres_0); u35_0.resize (nres_0);
  u36_0.resize  (nres_0); u37_0.resize  (nres_0); u38_0.resize (nres_0); u39_0.resize (nres_0); u40_0.resize (nres_0);
    
  for (int i = 0; i < nres_0; i++)
    { 
      mres_0(i) = mres_1(i);
      ntor_0(i) = ntor_1(i);
    }
 
  double weight1 = 1.;

  taua_0 = taua_1;
  P0_0   = P0_1;

  for (int j = 0; j < nres_0; j++)
    {
      u1_0 (j) = u1_1 (j); 
      u2_0 (j) = u2_1 (j); 
      u3_0 (j) = u3_1 (j); 
      u4_0 (j) = u4_1 (j); 
      u5_0 (j) = u5_1 (j); 
      u6_0 (j) = u6_1 (j); 
      u7_0 (j) = u7_1 (j); 
      u8_0 (j) = u8_1 (j); 
      u9_0 (j) = u9_1 (j); 
      u10_0(j) = u10_1(j); 
      u11_0(j) = u11_1(j); 
      u12_0(j) = u12_1(j); 
      u13_0(j) = u13_1(j); 
      u14_0(j) = u14_1(j); 
      u15_0(j) = u15_1(j); 
      u16_0(j) = u16_1(j); 
      u17_0(j) = u17_1(j); 
      u18_0(j) = u18_1(j); 
      u19_0(j) = u19_1(j);
      u20_0(j) = u20_1(j); 
      u21_0(j) = u21_1(j);
      u22_0(j) = u22_1(j); 
      u23_0(j) = u23_1(j); 
      u24_0(j) = u24_1(j); 
      u25_0(j) = u25_1(j);
      u26_0(j) = u26_1(j); 
      u27_0(j) = u27_1(j);
      u28_0(j) = u28_1(j); 
      u29_0(j) = u29_1(j);
      u30_0(j) = u30_1(j); 
      u31_0(j) = u31_1(j);
      u32_0(j) = u32_1(j); 
      u33_0(j) = u33_1(j); 
      u34_0(j) = u34_1(j); 
      u35_0(j) = u35_1(j);
      u36_0(j) = u36_1(j); 
      u37_0(j) = u37_1(j);
      u38_0(j) = u38_1(j);
      u39_0(j) = u39_1(j);
      u40_0(j) = u40_1(j); 
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e %16.9e\n", nres_0, taua_0, P0_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j),
	     u1_0(j),  u2_0(j),  u3_0(j),  u4_0(j),  u5_0(j),
	     u6_0(j),  u7_0(j),  u8_0(j),  u9_0(j),  u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j),
	     u16_0(j), u17_0(j), u18_0(j), u19_0(j), u20_0(j),
	     u21_0(j), u22_0(j), u23_0(j), u24_0(j), u25_0(j),
	     u26_0(j), u27_0(j), u28_0(j), u29_0(j), u30_0(j),
	     u31_0(j), u32_0(j), u33_0(j), u34_0(j), u35_0(j),
	     u36_0(j), u37_0(j), u38_0(j), u39_0(j), u40_0(j));
  
  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "%s %11.4e\n", nFile1, weight1);
  fclose (monitor);
}

void Phase::nFileInterpolateQuadratic (char* nFile1, double time1, char* nFile2, double time2, char* nFile, double time)
{
  int             nres_1;
  double          taua_1, P0_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1,  u2_1,  u3_1,  u4_1,  u5_1,  u6_1,  u7_1,  u8_1,  u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;
  Array<double,1> u20_1, u21_1, u22_1, u23_1, u24_1, u25_1, u26_1, u27_1, u28_1;
  Array<double,1> u29_1, u30_1, u31_1, u32_1, u33_1, u34_1, u35_1, u36_1, u37_1;
  Array<double,1> u38_1, u39_1, u40_1;

  int             nres_2;
  double          taua_2, P0_2;
  Array<int, 1>   mres_2, ntor_2;
  Array<double,1> u1_2,  u2_2,  u3_2,  u4_2,  u5_2,  u6_2,  u7_2,  u8_2,  u9_2, u10_2;
  Array<double,1> u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2, u18_2, u19_2;
  Array<double,1> u20_2, u21_2, u22_2, u23_2, u24_2, u25_2, u26_2, u27_2, u28_2;
  Array<double,1> u29_2, u30_2, u31_2, u32_2, u33_2, u34_2, u35_2, u36_2, u37_2;
  Array<double,1> u38_2, u39_2, u40_2;
       
  int             nres_0;
  double          taua_0, P0_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0,  u2_0,  u3_0,  u4_0,  u5_0,  u6_0,  u7_0,  u8_0,  u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
  Array<double,1> u20_0, u21_0, u22_0, u23_0, u24_0, u25_0, u26_0, u27_0, u28_0;
  Array<double,1> u29_0, u30_0, u31_0, u32_0, u33_0, u34_0, u35_0, u36_0, u37_0;
  Array<double,1> u38_0, u39_0, u40_0;
  
  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf %lf", &nres_1, &taua_1, &P0_1) != 3)
    {
      printf ("PHASE::nFileInterpolateQuadratic: Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1); ntor_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize   (nres_1); u3_1.resize  (nres_1); u4_1.resize  (nres_1); u5_1.resize  (nres_1);
  u6_1.resize   (nres_1); u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1); u10_1.resize (nres_1);
  u11_1.resize  (nres_1); u12_1.resize  (nres_1); u13_1.resize (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize  (nres_1); u18_1.resize (nres_1); u19_1.resize (nres_1); u20_1.resize (nres_1);
  u21_1.resize  (nres_1); u22_1.resize  (nres_1); u23_1.resize (nres_1); u24_1.resize (nres_1); u25_1.resize (nres_1);
  u26_1.resize  (nres_1); u27_1.resize  (nres_1); u28_1.resize (nres_1); u29_1.resize (nres_1); u30_1.resize (nres_1);
  u31_1.resize  (nres_1); u32_1.resize  (nres_1); u33_1.resize (nres_1); u34_1.resize (nres_1); u35_1.resize (nres_1);
  u36_1.resize  (nres_1); u37_1.resize  (nres_1); u38_1.resize (nres_1); u39_1.resize (nres_1); u40_1.resize (nres_1);
  
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j), &u20_1(j),
		  &u21_1(j), &u22_1(j), &u23_1(j), &u25_1(j), &u25_1(j),
		  &u26_1(j), &u27_1(j), &u28_1(j), &u29_1(j), &u30_1(j),
		  &u31_1(j), &u32_1(j), &u33_1(j), &u34_1(j), &u35_1(j),
		  &u36_1(j), &u37_1(j), &u38_1(j), &u39_1(j), &u40_1(j)) != 42)
	{
	  printf ("PHASE::nFileInterpolateQuadratic: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // .................
  // Read second nFile
  // .................
  file = OpenFiler (nFile2);
  
  if (fscanf (file, "%d %lf %lf", &nres_2, &taua_2, &P0_2) != 3)
    {
      printf ("PHASE::nFileInterpolateQuadratic: Error reading nFile_2 (1)\n");
      exit (1);
    }
  
  mres_2.resize (nres_2); ntor_2.resize (nres_2);
  u1_2.resize   (nres_2); u2_2.resize   (nres_2); u3_2.resize  (nres_2); u4_2.resize  (nres_2); u5_2.resize  (nres_2);
  u6_2.resize   (nres_2); u7_2.resize   (nres_2); u8_2.resize  (nres_2); u9_2.resize  (nres_2); u10_2.resize (nres_2);
  u11_2.resize  (nres_2); u12_2.resize  (nres_2); u13_2.resize (nres_2); u14_2.resize (nres_2); u15_2.resize (nres_2);
  u16_2.resize  (nres_2); u17_2.resize  (nres_2); u18_2.resize (nres_2); u19_2.resize (nres_2); u20_2.resize (nres_2);
  u21_2.resize  (nres_2); u22_2.resize  (nres_2); u23_2.resize (nres_2); u24_2.resize (nres_2); u25_2.resize (nres_2);
  u26_2.resize  (nres_2); u27_2.resize  (nres_2); u28_2.resize (nres_2); u29_2.resize (nres_2); u30_2.resize (nres_2);
  u31_2.resize  (nres_2); u32_2.resize  (nres_2); u33_2.resize (nres_2); u34_2.resize (nres_2); u35_2.resize (nres_2);
  u36_2.resize  (nres_2); u37_2.resize  (nres_2); u38_2.resize (nres_2); u39_2.resize (nres_2); u40_2.resize (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &ntor_2(j),
		  &u1_2 (j), &u2_2 (j), &u3_2 (j), &u4_2 (j), &u5_2 (j),
		  &u6_2 (j), &u7_2 (j), &u8_2 (j), &u9_2 (j), &u10_2(j),
		  &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j),
		  &u16_2(j), &u17_2(j), &u18_2(j), &u19_2(j), &u20_2(j),
		  &u21_2(j), &u22_2(j), &u23_2(j), &u25_2(j), &u25_2(j),
		  &u26_2(j), &u27_2(j), &u28_2(j), &u29_2(j), &u30_2(j),
		  &u31_2(j), &u32_2(j), &u33_2(j), &u34_2(j), &u35_2(j),
		  &u36_2(j), &u37_2(j), &u38_2(j), &u39_2(j), &u40_2(j)) != 42)
	{
	  printf ("PHASE::nFileInterpolateQuadratic: Error reading nFile_2 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ......................
  // Interpolate nFile data
  // ......................
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;

  mres_0.resize (nres_0); ntor_0.resize (nres_0);
  u1_0.resize   (nres_0); u2_0.resize   (nres_0); u3_0.resize  (nres_0); u4_0.resize  (nres_0); u5_0.resize  (nres_0);
  u6_0.resize   (nres_0); u7_0.resize   (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0); u10_0.resize (nres_0);
  u11_0.resize  (nres_0); u12_0.resize  (nres_0); u13_0.resize (nres_0); u19_0.resize (nres_0); u20_0.resize (nres_0);
  u21_0.resize  (nres_0); u22_0.resize  (nres_0); u23_0.resize (nres_0); u24_0.resize (nres_0); u25_0.resize (nres_0);
  u26_0.resize  (nres_0); u27_0.resize  (nres_0); u28_0.resize (nres_0); u29_0.resize (nres_0); u30_0.resize (nres_0);
  u31_0.resize  (nres_0); u32_0.resize  (nres_0); u33_0.resize (nres_0); u34_0.resize (nres_0); u35_0.resize (nres_0);
  u36_0.resize  (nres_0); u37_0.resize  (nres_0); u38_0.resize (nres_0); u39_0.resize (nres_0); u40_0.resize (nres_0);
  
  if (nres_0 == nres_1)
    {
      for (int i = 0; i < nres_0; i++)
        { 
	  mres_0(i) = mres_1(i);
	  ntor_0(i) = ntor_1(i);
	}
    }
  else 
    {
      for (int i = 0; i < nres_0; i++)
	{
	  mres_0(i) = mres_2(i);
	  ntor_0(i) = ntor_2(i);
	}
    }
  
  double weight1 = (time - time2)  /(time1 - time2);
  double weight2 = (time - time1)  /(time2 - time1);

  taua_0 = weight1 * taua_1 + weight2 * taua_2;
  P0_0   = weight1 * P0_1   + weight2 * P0_2;

  for (int j = 0; j < nres_0; j++)
    {
      u1_0 (j) = weight1 * u1_1 (j) + weight2 * u1_2 (j);
      u2_0 (j) = weight1 * u2_1 (j) + weight2 * u2_2 (j);
      u3_0 (j) = weight1 * u3_1 (j) + weight2 * u3_2 (j);
      u4_0 (j) = weight1 * u4_1 (j) + weight2 * u4_2 (j);
      u5_0 (j) = weight1 * u5_1 (j) + weight2 * u5_2 (j);
      u6_0 (j) = weight1 * u6_1 (j) + weight2 * u6_2 (j);
      u7_0 (j) = weight1 * u7_1 (j) + weight2 * u7_2 (j);
      u8_0 (j) = weight1 * u8_1 (j) + weight2 * u8_2 (j);
      u9_0 (j) = weight1 * u9_1 (j) + weight2 * u9_2 (j);
      u10_0(j) = weight1 * u10_1(j) + weight2 * u10_2(j);
      u11_0(j) = weight1 * u11_1(j) + weight2 * u11_2(j);
      u12_0(j) = weight1 * u12_1(j) + weight2 * u12_2(j);
      u13_0(j) = weight1 * u13_1(j) + weight2 * u13_2(j);
      u14_0(j) = weight1 * u14_1(j) + weight2 * u14_2(j);
      u15_0(j) = weight1 * u15_1(j) + weight2 * u15_2(j);
      u16_0(j) = weight1 * u16_1(j) + weight2 * u16_2(j);
      u17_0(j) = weight1 * u17_1(j) + weight2 * u17_2(j);
      u18_0(j) = weight1 * u18_1(j) + weight2 * u18_2(j);
      u19_0(j) = weight1 * u19_1(j) + weight2 * u19_2(j);
      u20_0(j) = weight1 * u20_1(j) + weight2 * u20_2(j);
      u21_0(j) = weight1 * u21_1(j) + weight2 * u21_2(j);
      u22_0(j) = weight1 * u22_1(j) + weight2 * u22_2(j);
      u23_0(j) = weight1 * u23_1(j) + weight2 * u23_2(j);
      u24_0(j) = weight1 * u24_1(j) + weight2 * u24_2(j);
      u25_0(j) = weight1 * u25_1(j) + weight2 * u25_2(j);
      u26_0(j) = weight1 * u26_1(j) + weight2 * u26_2(j);
      u27_0(j) = weight1 * u27_1(j) + weight2 * u27_2(j);
      u28_0(j) = weight1 * u28_1(j) + weight2 * u28_2(j);
      u29_0(j) = weight1 * u29_1(j) + weight2 * u29_2(j);
      u30_0(j) = weight1 * u30_1(j) + weight2 * u30_2(j);
      u31_0(j) = weight1 * u31_1(j) + weight2 * u31_2(j);
      u32_0(j) = weight1 * u32_1(j) + weight2 * u32_2(j);
      u33_0(j) = weight1 * u33_1(j) + weight2 * u33_2(j);
      u34_0(j) = weight1 * u34_1(j) + weight2 * u34_2(j);
      u35_0(j) = weight1 * u35_1(j) + weight2 * u35_2(j);
      u36_0(j) = weight1 * u36_1(j) + weight2 * u36_2(j);
      u37_0(j) = weight1 * u37_1(j) + weight2 * u37_2(j);
      u38_0(j) = weight1 * u38_1(j) + weight2 * u38_2(j);
      u39_0(j) = weight1 * u39_1(j) + weight2 * u39_2(j);
      u40_0(j) = weight1 * u40_1(j) + weight2 * u40_2(j);
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e %16.9e\n", nres_0, taua_0, P0_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j),
	     u1_0(j),  u2_0(j),  u3_0(j),  u4_0(j),  u5_0(j),
	     u6_0(j),  u7_0(j),  u8_0(j),  u9_0(j),  u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j),
	     u16_0(j), u17_0(j), u18_0(j), u19_0(j), u20_0(j),
	     u21_0(j), u22_0(j), u23_0(j), u24_0(j), u25_0(j),
	     u26_0(j), u27_0(j), u28_0(j), u29_0(j), u30_0(j),
	     u31_0(j), u32_0(j), u33_0(j), u34_0(j), u35_0(j),
	     u36_0(j), u37_0(j), u38_0(j), u39_0(j), u40_0(j));
      
  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);
  printf ("%s %11.4e\n", nFile2, weight2);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "%s %11.4e\n", nFile1, weight1);
  fprintf (monitor, "%s %11.4e\n", nFile2, weight2);
  fclose (monitor);
}
  
void Phase::nFileInterpolateCubic (char* nFile1, double time1, char* nFile2, double time2, char* nFile3, double time3, char* nFile, double time)
{
  int             nres_1;
  double          taua_1, P0_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1,  u2_1,  u3_1,  u4_1,  u5_1,  u6_1,  u7_1,  u8_1,  u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;
  Array<double,1> u20_1, u21_1, u22_1, u23_1, u24_1, u25_1, u26_1, u27_1, u28_1;
  Array<double,1> u29_1, u30_1, u31_1, u32_1, u33_1, u34_1, u35_1, u36_1, u37_1;
  Array<double,1> u38_1, u39_1, u40_1;

  int             nres_2;
  double          taua_2, P0_2;
  Array<int, 1>   mres_2, ntor_2;
  Array<double,1> u1_2,  u2_2,  u3_2,  u4_2,  u5_2,  u6_2,  u7_2,  u8_2,  u9_2, u10_2;
  Array<double,1> u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2, u18_2, u19_2;
  Array<double,1> u20_2, u21_2, u22_2, u23_2, u24_2, u25_2, u26_2, u27_2, u28_2;
  Array<double,1> u29_2, u30_2, u31_2, u32_2, u33_2, u34_2, u35_2, u36_2, u37_2;
  Array<double,1> u38_2, u39_2, u40_2;

  int             nres_3;
  double          taua_3, P0_3;
  Array<int, 1>   mres_3, ntor_3;
  Array<double,1> u1_3,  u2_3,  u3_3,  u4_3,  u5_3,  u6_3,  u7_3,  u8_3,  u9_3, u10_3;
  Array<double,1> u11_3, u12_3, u13_3, u14_3, u15_3, u16_3, u17_3, u18_3, u19_3;
  Array<double,1> u20_3, u21_3, u22_3, u23_3, u24_3, u25_3, u26_3, u27_3, u28_3;
  Array<double,1> u29_3, u30_3, u31_3, u32_3, u33_3, u34_3, u35_3, u36_3, u37_3;
  Array<double,1> u38_3, u39_3, u40_3;
         
  int             nres_0;
  double          taua_0, P0_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0,  u2_0,  u3_0,  u4_0,  u5_0,  u6_0,  u7_0,  u8_0,  u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
  Array<double,1> u20_0, u21_0, u22_0, u23_0, u24_0, u25_0, u26_0, u27_0, u28_0;
  Array<double,1> u29_0, u30_0, u31_0, u32_0, u33_0, u34_0, u35_0, u36_0, u37_0;
  Array<double,1> u38_0, u39_0, u40_0;

  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf %lf", &nres_1, &taua_1, &P0_1) != 3)
    {
      printf ("PHASE::nFileInterpolateCubic: Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1); ntor_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize   (nres_1); u3_1.resize  (nres_1); u4_1.resize  (nres_1); u5_1.resize  (nres_1);
  u6_1.resize   (nres_1); u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1); u10_1.resize (nres_1);
  u11_1.resize  (nres_1); u12_1.resize  (nres_1); u13_1.resize (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize  (nres_1); u18_1.resize (nres_1); u19_1.resize (nres_1); u20_1.resize (nres_1);
  u21_1.resize  (nres_1); u22_1.resize  (nres_1); u23_1.resize (nres_1); u24_1.resize (nres_1); u25_1.resize (nres_1);
  u26_1.resize  (nres_1); u27_1.resize  (nres_1); u28_1.resize (nres_1); u29_1.resize (nres_1); u30_1.resize (nres_1);
  u31_1.resize  (nres_1); u32_1.resize  (nres_1); u33_1.resize (nres_1); u34_1.resize (nres_1); u35_1.resize (nres_1);
  u36_1.resize  (nres_1); u37_1.resize  (nres_1); u38_1.resize (nres_1); u39_1.resize (nres_1); u40_1.resize (nres_1);
    
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j), &u20_1(j),
		  &u21_1(j), &u22_1(j), &u23_1(j), &u25_1(j), &u25_1(j),
		  &u26_1(j), &u27_1(j), &u28_1(j), &u29_1(j), &u30_1(j),
		  &u31_1(j), &u32_1(j), &u33_1(j), &u34_1(j), &u35_1(j),
		  &u36_1(j), &u37_1(j), &u38_1(j), &u39_1(j), &u40_1(j)) != 42)
	{
	  printf ("PHASE: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // .................
  // Read second nFile
  // .................
  file = OpenFiler (nFile2);
  
  if (fscanf (file, "%d %lf %lf", &nres_2, &taua_2, &P0_2) != 3)
    {
      printf ("PHASE::nFileInterpolateCubic: Error reading nFile_2 (1)\n");
      exit (1);
    }
  
  mres_2.resize (nres_2); ntor_2.resize (nres_2);
  u1_2.resize   (nres_2); u2_2.resize   (nres_2); u3_2.resize  (nres_2); u4_2.resize  (nres_2); u5_2.resize  (nres_2);
  u6_2.resize   (nres_2); u7_2.resize   (nres_2); u8_2.resize  (nres_2); u9_2.resize  (nres_2); u10_2.resize (nres_2);
  u11_2.resize  (nres_2); u12_2.resize  (nres_2); u13_2.resize (nres_2); u14_2.resize (nres_2); u15_2.resize (nres_2);
  u16_2.resize  (nres_2); u17_2.resize  (nres_2); u18_2.resize (nres_2); u19_2.resize (nres_2); u20_2.resize (nres_2);
  u21_2.resize  (nres_2); u22_2.resize  (nres_2); u23_2.resize (nres_2); u24_2.resize (nres_2); u25_2.resize (nres_2);
  u26_2.resize  (nres_2); u27_2.resize  (nres_2); u28_2.resize (nres_2); u29_2.resize (nres_2); u30_2.resize (nres_2);
  u31_2.resize  (nres_2); u32_2.resize  (nres_2); u33_2.resize (nres_2); u34_2.resize (nres_2); u35_2.resize (nres_2);
  u36_2.resize  (nres_2); u37_2.resize  (nres_2); u38_2.resize (nres_2); u39_2.resize (nres_2); u40_2.resize (nres_2);
   
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &ntor_2(j),
		  &u1_2 (j), &u2_2 (j), &u3_2 (j), &u4_2 (j), &u5_2 (j),
		  &u6_2 (j), &u7_2 (j), &u8_2 (j), &u9_2 (j), &u10_2(j),
		  &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j),
		  &u16_2(j), &u17_2(j), &u18_2(j), &u19_2(j), &u20_2(j),
		  &u21_2(j), &u22_2(j), &u23_2(j), &u25_2(j), &u25_2(j),
		  &u26_2(j), &u27_2(j), &u28_2(j), &u29_2(j), &u30_2(j),
		  &u31_2(j), &u32_2(j), &u33_2(j), &u34_2(j), &u35_2(j),
		  &u36_2(j), &u37_2(j), &u38_2(j), &u39_2(j), &u40_2(j)) != 42)
	{
	  printf ("PHASE::nFileInterpolateCubic: Error reading nFile_2 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ................
  // Read third nFile
  // ................
  file = OpenFiler (nFile3);

  if (fscanf (file, "%d %lf %lf", &nres_3, &taua_3, &P0_3) != 3)
    {
      printf ("PHASE::nFileInterpolateCubic: Error reading nFile_3 (1)\n");
      exit (1);
    }

  mres_3.resize (nres_3); ntor_3.resize (nres_3);
  u1_3.resize   (nres_3); u2_3.resize   (nres_3); u3_3.resize  (nres_3); u4_3.resize  (nres_3); u5_3.resize  (nres_3);
  u6_3.resize   (nres_3); u7_3.resize   (nres_3); u8_3.resize  (nres_3); u9_3.resize  (nres_3); u10_3.resize (nres_3);
  u11_3.resize  (nres_3); u12_3.resize  (nres_3); u13_3.resize (nres_3); u14_3.resize (nres_3); u15_3.resize (nres_3);
  u16_3.resize  (nres_3); u17_3.resize  (nres_3); u18_3.resize (nres_3); u19_3.resize (nres_3); u20_3.resize (nres_3);
  u21_3.resize  (nres_3); u22_3.resize  (nres_3); u23_3.resize (nres_3); u24_3.resize (nres_3); u25_3.resize (nres_3);
  u26_3.resize  (nres_3); u27_3.resize  (nres_3); u28_3.resize (nres_3); u29_3.resize (nres_3); u30_3.resize (nres_3);
  u31_3.resize  (nres_3); u32_3.resize  (nres_3); u33_3.resize (nres_3); u34_3.resize (nres_3); u35_3.resize (nres_3);
  u36_3.resize  (nres_3); u37_3.resize  (nres_3); u38_3.resize (nres_3); u39_3.resize (nres_3); u40_3.resize (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &ntor_3(j),
		  &u1_3 (j), &u2_3 (j), &u3_3 (j), &u4_3 (j), &u5_3 (j),
		  &u6_3 (j), &u7_3 (j), &u8_3 (j), &u9_3 (j), &u10_3(j),
		  &u11_3(j), &u12_3(j), &u13_3(j), &u14_3(j), &u15_3(j),
		  &u16_3(j), &u17_3(j), &u18_3(j), &u19_3(j), &u20_3(j),
		  &u21_3(j), &u22_3(j), &u23_3(j), &u25_3(j), &u25_3(j),
		  &u26_3(j), &u27_3(j), &u28_3(j), &u29_3(j), &u30_3(j),
		  &u31_3(j), &u32_3(j), &u33_3(j), &u34_3(j), &u35_3(j),
		  &u36_3(j), &u37_3(j), &u38_3(j), &u39_3(j), &u40_3(j)) != 42)
 	{
	  printf ("PHASE::nFileInterpolateCubic: Error reading nFile_3 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ......................
  // Interpolate nFile data
  // ......................
  if (nres_1 == nres_2 && nres_2 != nres_3)
    {
      nFileInterpolateQuadratic (nFile1, time1, nFile2, time2, nFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2)
    {
      nFileInterpolateQuadratic (nFile2, time2, nFile3, time3, nFile, time);
      return;
    }
  
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
  if (nres_3 < nres_0)
    nres_0 = nres_3;
  
  mres_0.resize (nres_0); ntor_0.resize (nres_0);
  u1_0.resize   (nres_0); u2_0.resize   (nres_0); u3_0.resize  (nres_0); u4_0.resize  (nres_0); u5_0.resize  (nres_0);
  u6_0.resize   (nres_0); u7_0.resize   (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0); u10_0.resize (nres_0);
  u11_0.resize  (nres_0); u12_0.resize  (nres_0); u13_0.resize (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize  (nres_0); u17_0.resize  (nres_0); u18_0.resize (nres_0); u19_0.resize (nres_0); u20_0.resize (nres_0);
  u21_0.resize  (nres_0); u22_0.resize  (nres_0); u23_0.resize (nres_0); u24_0.resize (nres_0); u25_0.resize (nres_0);
  u26_0.resize  (nres_0); u27_0.resize  (nres_0); u28_0.resize (nres_0); u29_0.resize (nres_0); u30_0.resize (nres_0);
  u31_0.resize  (nres_0); u32_0.resize  (nres_0); u33_0.resize (nres_0); u34_0.resize (nres_0); u35_0.resize (nres_0);
  u36_0.resize  (nres_0); u37_0.resize  (nres_0); u38_0.resize (nres_0); u39_0.resize (nres_0); u40_0.resize (nres_0);
  
  if (nres_0 == nres_1)
    {
      for (int i = 0; i < nres_0; i++)
        { 
	  mres_0(i) = mres_1(i);
	  ntor_0(i) = ntor_1(i);
	}
    }
  else if (nres_0 == nres_2)
    {
      for (int i = 0; i < nres_0; i++)
	{
	  mres_0(i) = mres_2(i);
	  ntor_0(i) = ntor_2(i);
	}
    }
  else
    {
      for (int i = 0; i < nres_0; i++)
	{
	  mres_0(i) = mres_3(i);
	  ntor_0(i) = ntor_3(i);
	}
    }
  
  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);

  taua_0 = weight1 * taua_1 + weight2 * taua_2 + weight3 * taua_3;
  P0_0   = weight1 * P0_1   + weight2 * P0_2   + weight3 * P0_3;

  for (int j = 0; j < nres_0; j++)
    {
      u1_0 (j) = weight1 * u1_1 (j) + weight2 * u1_2 (j) + weight3 * u1_3 (j);
      u2_0 (j) = weight1 * u2_1 (j) + weight2 * u2_2 (j) + weight3 * u2_3 (j);
      u3_0 (j) = weight1 * u3_1 (j) + weight2 * u3_2 (j) + weight3 * u3_3 (j);
      u4_0 (j) = weight1 * u4_1 (j) + weight2 * u4_2 (j) + weight3 * u4_3 (j);
      u5_0 (j) = weight1 * u5_1 (j) + weight2 * u5_2 (j) + weight3 * u5_3 (j);
      u6_0 (j) = weight1 * u6_1 (j) + weight2 * u6_2 (j) + weight3 * u6_3 (j);
      u7_0 (j) = weight1 * u7_1 (j) + weight2 * u7_2 (j) + weight3 * u7_3 (j);
      u8_0 (j) = weight1 * u8_1 (j) + weight2 * u8_2 (j) + weight3 * u8_3 (j);
      u9_0 (j) = weight1 * u9_1 (j) + weight2 * u9_2 (j) + weight3 * u9_3 (j);
      u10_0(j) = weight1 * u10_1(j) + weight2 * u10_2(j) + weight3 * u10_3(j);
      u11_0(j) = weight1 * u11_1(j) + weight2 * u11_2(j) + weight3 * u11_3(j);
      u12_0(j) = weight1 * u12_1(j) + weight2 * u12_2(j) + weight3 * u12_3(j);
      u13_0(j) = weight1 * u13_1(j) + weight2 * u13_2(j) + weight3 * u13_3(j);
      u14_0(j) = weight1 * u14_1(j) + weight2 * u14_2(j) + weight3 * u14_3(j);
      u15_0(j) = weight1 * u15_1(j) + weight2 * u15_2(j) + weight3 * u15_3(j);
      u16_0(j) = weight1 * u16_1(j) + weight2 * u16_2(j) + weight3 * u16_3(j);
      u17_0(j) = weight1 * u17_1(j) + weight2 * u17_2(j) + weight3 * u17_3(j);
      u18_0(j) = weight1 * u18_1(j) + weight2 * u18_2(j) + weight3 * u18_3(j);
      u19_0(j) = weight1 * u19_1(j) + weight2 * u19_2(j) + weight3 * u19_3(j);
      u20_0(j) = weight1 * u20_1(j) + weight2 * u20_2(j) + weight3 * u20_3(j);
      u21_0(j) = weight1 * u21_1(j) + weight2 * u21_2(j) + weight3 * u21_3(j);
      u22_0(j) = weight1 * u22_1(j) + weight2 * u22_2(j) + weight3 * u22_3(j);
      u23_0(j) = weight1 * u23_1(j) + weight2 * u23_2(j) + weight3 * u23_3(j);
      u24_0(j) = weight1 * u24_1(j) + weight2 * u24_2(j) + weight3 * u24_3(j);
      u25_0(j) = weight1 * u25_1(j) + weight2 * u25_2(j) + weight3 * u25_3(j);
      u26_0(j) = weight1 * u26_1(j) + weight2 * u26_2(j) + weight3 * u26_3(j);
      u27_0(j) = weight1 * u27_1(j) + weight2 * u27_2(j) + weight3 * u27_3(j);
      u28_0(j) = weight1 * u28_1(j) + weight2 * u28_2(j) + weight3 * u28_3(j);
      u29_0(j) = weight1 * u29_1(j) + weight2 * u29_2(j) + weight3 * u29_3(j);
      u30_0(j) = weight1 * u30_1(j) + weight2 * u30_2(j) + weight3 * u30_3(j);
      u31_0(j) = weight1 * u31_1(j) + weight2 * u31_2(j) + weight3 * u31_3(j);
      u32_0(j) = weight1 * u32_1(j) + weight2 * u32_2(j) + weight3 * u32_3(j);
      u33_0(j) = weight1 * u33_1(j) + weight2 * u33_2(j) + weight3 * u33_3(j);
      u34_0(j) = weight1 * u34_1(j) + weight2 * u34_2(j) + weight3 * u34_3(j);
      u35_0(j) = weight1 * u35_1(j) + weight2 * u35_2(j) + weight3 * u35_3(j);
      u36_0(j) = weight1 * u36_1(j) + weight2 * u36_2(j) + weight3 * u36_3(j);
      u37_0(j) = weight1 * u37_1(j) + weight2 * u37_2(j) + weight3 * u37_3(j);
      u38_0(j) = weight1 * u38_1(j) + weight2 * u38_2(j) + weight3 * u38_3(j);
      u39_0(j) = weight1 * u39_1(j) + weight2 * u39_2(j) + weight3 * u39_3(j);
      u40_0(j) = weight1 * u40_1(j) + weight2 * u40_2(j) + weight3 * u40_3(j);
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e %16.9e\n", nres_0, taua_0, P0_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j),
	     u1_0(j),  u2_0(j),  u3_0(j),  u4_0(j),  u5_0(j),
	     u6_0(j),  u7_0(j),  u8_0(j),  u9_0(j),  u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j),
	     u16_0(j), u17_0(j), u18_0(j), u19_0(j), u20_0(j),
	     u21_0(j), u22_0(j), u23_0(j), u24_0(j), u25_0(j),
	     u26_0(j), u27_0(j), u28_0(j), u29_0(j), u30_0(j),
	     u31_0(j), u32_0(j), u33_0(j), u34_0(j), u35_0(j),
	     u36_0(j), u37_0(j), u38_0(j), u39_0(j), u40_0(j));
  
  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);
  printf ("%s %11.4e\n", nFile2, weight2);
  printf ("%s %11.4e\n", nFile3, weight3);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "%s %11.4e\n", nFile1, weight1);
  fprintf (monitor, "%s %11.4e\n", nFile2, weight2);
  fprintf (monitor, "%s %11.4e\n", nFile3, weight3);
  fclose (monitor);
}

void Phase::nFileInterpolateQuartic (char* nFile1, double time1, char* nFile2, double time2, char* nFile3, double time3,
				     char* nFile4, double time4, char* nFile, double time)
{
  int             nres_1;
  double          taua_1, P0_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1,  u2_1,  u3_1,  u4_1,  u5_1,  u6_1,  u7_1,  u8_1,  u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;
  Array<double,1> u20_1, u21_1, u22_1, u23_1, u24_1, u25_1, u26_1, u27_1, u28_1;
  Array<double,1> u29_1, u30_1, u31_1, u32_1, u33_1, u34_1, u35_1, u36_1, u37_1;
  Array<double,1> u38_1, u39_1, u40_1;

  int             nres_2;
  double          taua_2, P0_2;
  Array<int, 1>   mres_2, ntor_2;
  Array<double,1> u1_2,  u2_2,  u3_2,  u4_2,  u5_2,  u6_2,  u7_2,  u8_2,  u9_2, u10_2;
  Array<double,1> u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2, u18_2, u19_2;
  Array<double,1> u20_2, u21_2, u22_2, u23_2, u24_2, u25_2, u26_2, u27_2, u28_2;
  Array<double,1> u29_2, u30_2, u31_2, u32_2, u33_2, u34_2, u35_2, u36_2, u37_2;
  Array<double,1> u38_2, u39_2, u40_2;

  int             nres_3;
  double          taua_3, P0_3;
  Array<int, 1>   mres_3, ntor_3;
  Array<double,1> u1_3,  u2_3,  u3_3,  u4_3,  u5_3,  u6_3,  u7_3,  u8_3,  u9_3, u10_3;
  Array<double,1> u11_3, u12_3, u13_3, u14_3, u15_3, u16_3, u17_3, u18_3, u19_3;
  Array<double,1> u20_3, u21_3, u22_3, u23_3, u24_3, u25_3, u26_3, u27_3, u28_3;
  Array<double,1> u29_3, u30_3, u31_3, u32_3, u33_3, u34_3, u35_3, u36_3, u37_3;
  Array<double,1> u38_3, u39_3, u40_3;

  int             nres_4;
  double          taua_4, P0_4;
  Array<int, 1>   mres_4, ntor_4;
  Array<double,1> u1_4,  u2_4,  u3_4,  u4_4,  u5_4,  u6_4,  u7_4,  u8_4,  u9_4, u10_4;
  Array<double,1> u11_4, u12_4, u13_4, u14_4, u15_4, u16_4, u17_4, u18_4, u19_4;
  Array<double,1> u20_4, u21_4, u22_4, u23_4, u24_4, u25_4, u26_4, u27_4, u28_4;
  Array<double,1> u29_4, u30_4, u31_4, u32_4, u33_4, u34_4, u35_4, u36_4, u37_4;
  Array<double,1> u38_4, u39_4, u40_4;
         
  int             nres_0;
  double          taua_0, P0_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0,  u2_0,  u3_0,  u4_0,  u5_0,  u6_0,  u7_0,  u8_0,  u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
  Array<double,1> u20_0, u21_0, u22_0, u23_0, u24_0, u25_0, u26_0, u27_0, u28_0;
  Array<double,1> u29_0, u30_0, u31_0, u32_0, u33_0, u34_0, u35_0, u36_0, u37_0;
  Array<double,1> u38_0, u39_0, u40_0;

  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf %lf", &nres_1, &taua_1, &P0_1) != 3)
    {
      printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1); ntor_1.resize (nres_1);
  u1_1.resize   (nres_1); u2_1.resize   (nres_1); u3_1.resize  (nres_1); u4_1.resize  (nres_1); u5_1.resize  (nres_1);
  u6_1.resize   (nres_1); u7_1.resize   (nres_1); u8_1.resize  (nres_1); u9_1.resize  (nres_1); u10_1.resize (nres_1);
  u11_1.resize  (nres_1); u12_1.resize  (nres_1); u13_1.resize (nres_1); u14_1.resize (nres_1); u15_1.resize (nres_1);
  u16_1.resize  (nres_1); u17_1.resize  (nres_1); u18_1.resize (nres_1); u19_1.resize (nres_1); u20_1.resize (nres_1);
  u21_1.resize  (nres_1); u22_1.resize  (nres_1); u23_1.resize (nres_1); u24_1.resize (nres_1); u25_1.resize (nres_1);
  u26_1.resize  (nres_1); u27_1.resize  (nres_1); u28_1.resize (nres_1); u29_1.resize (nres_1); u30_1.resize (nres_1);
  u31_1.resize  (nres_1); u32_1.resize  (nres_1); u33_1.resize (nres_1); u34_1.resize (nres_1); u35_1.resize (nres_1);
  u36_1.resize  (nres_1); u37_1.resize  (nres_1); u38_1.resize (nres_1); u39_1.resize (nres_1); u40_1.resize (nres_1);
    
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j), &u20_1(j),
		  &u21_1(j), &u22_1(j), &u23_1(j), &u25_1(j), &u25_1(j),
		  &u26_1(j), &u27_1(j), &u28_1(j), &u29_1(j), &u30_1(j),
		  &u31_1(j), &u32_1(j), &u33_1(j), &u34_1(j), &u35_1(j),
		  &u36_1(j), &u37_1(j), &u38_1(j), &u39_1(j), &u40_1(j)) != 42)
      	{
	  printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // .................
  // Read second nFile
  // .................
  file = OpenFiler (nFile2);
  
  if (fscanf (file, "%d %lf %lf", &nres_2, &taua_2, &P0_2) != 3)
    {
      printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_2 (1)\n");
      exit (1);
    }
  
  mres_2.resize (nres_2); ntor_2.resize (nres_2);
  u1_2.resize   (nres_2); u2_2.resize   (nres_2); u3_2.resize  (nres_2); u4_2.resize   (nres_2); u5_2.resize  (nres_2);
  u6_2.resize   (nres_2); u7_2.resize   (nres_2); u8_2.resize  (nres_2); u9_2.resize   (nres_2); u10_2.resize (nres_2);
  u11_2.resize  (nres_2); u12_2.resize  (nres_2); u13_2.resize (nres_2); u14_2.resize  (nres_2); u15_2.resize (nres_2);
  u16_2.resize  (nres_2); u17_2.resize  (nres_2); u18_2.resize (nres_2); u19_2.resize  (nres_2); u20_2.resize (nres_2);
  u21_2.resize  (nres_2); u22_2.resize  (nres_2); u23_2.resize (nres_2); u24_2.resize  (nres_2); u25_2.resize (nres_2);
  u26_2.resize  (nres_2); u27_2.resize  (nres_2); u28_2.resize (nres_2); u29_2.resize  (nres_2); u30_2.resize (nres_2);
  u31_2.resize  (nres_2); u32_2.resize  (nres_2); u33_2.resize (nres_2); u34_2.resize  (nres_2); u35_2.resize (nres_2);
  u36_2.resize  (nres_2); u37_2.resize  (nres_2); u38_2.resize (nres_2); u39_2.resize  (nres_2); u40_2.resize (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &ntor_2(j),
		  &u1_2 (j), &u2_2 (j), &u3_2 (j), &u4_2 (j), &u5_2 (j),
		  &u6_2 (j), &u7_2 (j), &u8_2 (j), &u9_2 (j), &u10_2(j),
		  &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j),
		  &u16_2(j), &u17_2(j), &u18_2(j), &u19_2(j), &u20_2(j),
		  &u21_2(j), &u22_2(j), &u23_2(j), &u25_2(j), &u25_2(j),
		  &u26_2(j), &u27_2(j), &u28_2(j), &u29_2(j), &u30_2(j),
		  &u31_2(j), &u32_2(j), &u33_2(j), &u34_2(j), &u35_2(j),
		  &u36_2(j), &u37_2(j), &u38_2(j), &u39_2(j), &u40_2(j)) != 42)
     	{
	  printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_2 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ................
  // Read third nFile
  // ................
  file = OpenFiler (nFile3);

  if (fscanf (file, "%d %lf %lf", &nres_3, &taua_3, &P0_3) != 3)
    {
      printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_3 (1)\n");
      exit (1);
    }

  mres_3.resize (nres_3); ntor_3.resize (nres_3);
  u1_3.resize   (nres_3); u2_3.resize   (nres_3); u3_3.resize  (nres_3); u4_3.resize  (nres_3); u5_3.resize  (nres_3);
  u6_3.resize   (nres_3); u7_3.resize   (nres_3); u8_3.resize  (nres_3); u9_3.resize  (nres_3); u10_3.resize (nres_3);
  u11_3.resize  (nres_3); u12_3.resize  (nres_3); u13_3.resize (nres_3); u14_3.resize (nres_3); u15_3.resize (nres_3);
  u16_3.resize  (nres_3); u17_3.resize  (nres_3); u18_3.resize (nres_3); u19_3.resize (nres_3); u20_3.resize (nres_3);
  u21_3.resize  (nres_3); u22_3.resize  (nres_3); u23_3.resize (nres_3); u24_3.resize (nres_3); u25_3.resize (nres_3);
  u26_3.resize  (nres_3); u27_3.resize  (nres_3); u28_3.resize (nres_3); u29_3.resize (nres_3); u30_3.resize (nres_3);
  u31_3.resize  (nres_3); u32_3.resize  (nres_3); u33_3.resize (nres_3); u34_3.resize (nres_3); u35_3.resize (nres_3);
  u36_3.resize  (nres_3); u37_3.resize  (nres_3); u38_3.resize (nres_3); u39_3.resize (nres_3); u40_3.resize (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &ntor_3(j),
		  &u1_3 (j), &u2_3 (j), &u3_3 (j), &u4_3 (j), &u5_3 (j),
		  &u6_3 (j), &u7_3 (j), &u8_3 (j), &u9_3 (j), &u10_3(j),
		  &u11_3(j), &u12_3(j), &u13_3(j), &u14_3(j), &u15_3(j),
		  &u16_3(j), &u17_3(j), &u18_3(j), &u19_3(j), &u20_3(j),
		  &u21_3(j), &u22_3(j), &u23_3(j), &u25_3(j), &u25_3(j),
		  &u26_3(j), &u27_3(j), &u28_3(j), &u29_3(j), &u30_3(j),
		  &u31_3(j), &u32_3(j), &u33_3(j), &u34_3(j), &u35_3(j),
		  &u36_3(j), &u37_3(j), &u38_3(j), &u39_3(j), &u40_3(j)) != 42)
	{
	  printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_3 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);

  // .................
  // Read fourth nFile
  // .................
  file = OpenFiler (nFile4);

  if (fscanf (file, "%d %lf %lf", &nres_4, &taua_4, &P0_4) != 3)
    {
      printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_4 (1)\n");
      exit (1);
    }

  mres_4.resize (nres_4); ntor_4.resize (nres_4);
  u1_4.resize   (nres_4); u2_4.resize   (nres_4); u3_4.resize  (nres_4); u4_4.resize  (nres_4); u5_4.resize  (nres_4);
  u6_4.resize   (nres_4); u7_4.resize   (nres_4); u8_4.resize  (nres_4); u9_4.resize  (nres_4); u10_4.resize (nres_4);
  u11_4.resize  (nres_4); u12_4.resize  (nres_4); u13_4.resize (nres_4); u14_4.resize (nres_4); u15_4.resize (nres_4);
  u16_4.resize  (nres_4); u17_4.resize  (nres_4); u18_4.resize (nres_4); u19_4.resize (nres_4); u20_4.resize (nres_4);
  u21_4.resize  (nres_4); u22_4.resize  (nres_4); u23_4.resize (nres_4); u24_4.resize (nres_4); u25_4.resize (nres_4);
  u26_4.resize  (nres_4); u27_4.resize  (nres_4); u28_4.resize (nres_4); u29_4.resize (nres_4); u30_4.resize (nres_4);
  u31_4.resize  (nres_4); u32_4.resize  (nres_4); u33_4.resize (nres_4); u34_4.resize (nres_4); u35_4.resize (nres_4);
  u36_4.resize  (nres_4); u37_4.resize  (nres_4); u38_4.resize (nres_4); u39_4.resize (nres_4); u40_4.resize (nres_4);
  
  for (int j = 0; j < nres_4; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_4(j), &ntor_4(j),
		  &u1_4 (j), &u2_4 (j), &u3_4 (j), &u4_4 (j), &u5_4 (j),
		  &u6_4 (j), &u7_4 (j), &u8_4 (j), &u9_4 (j), &u10_4(j),
		  &u11_4(j), &u12_4(j), &u13_4(j), &u14_4(j), &u15_4(j),
		  &u16_4(j), &u17_4(j), &u18_4(j), &u19_4(j), &u20_4(j),
		  &u21_4(j), &u22_4(j), &u23_4(j), &u25_4(j), &u25_4(j),
		  &u26_4(j), &u27_4(j), &u28_4(j), &u29_4(j), &u30_4(j),
		  &u31_4(j), &u32_4(j), &u33_4(j), &u34_4(j), &u35_4(j),
		  &u36_4(j), &u37_4(j), &u38_4(j), &u39_4(j), &u40_4(j)) != 42)
	{
	  printf ("PHASE::nFileInterpolateQuartic: Error reading nFile_4 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ......................
  // Interpolate nFile data
  // ......................
  if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 != nres_2)
    {
      nFileInterpolateQuadratic (nFile2, time2, nFile3, time3, nFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 == nres_2 && nres_4 != nres_2)
    {
      nFileInterpolateCubic (nFile1, time1, nFile2, time2, nFile3, time3, nFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 == nres_2)
    {
      nFileInterpolateCubic (nFile2, time2, nFile3, time3, nFile4, time4, nFile, time);
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
  
  mres_0.resize (nres_0); ntor_0.resize (nres_0);
  u1_0.resize   (nres_0); u2_0.resize   (nres_0); u3_0.resize  (nres_0); u4_0.resize  (nres_0); u5_0.resize  (nres_0);
  u6_0.resize   (nres_0); u7_0.resize   (nres_0); u8_0.resize  (nres_0); u9_0.resize  (nres_0); u10_0.resize (nres_0);
  u11_0.resize  (nres_0); u12_0.resize  (nres_0); u13_0.resize (nres_0); u14_0.resize (nres_0); u15_0.resize (nres_0);
  u16_0.resize  (nres_0); u17_0.resize  (nres_0); u18_0.resize (nres_0); u19_0.resize (nres_0); u20_0.resize (nres_0);
  u21_0.resize  (nres_0); u22_0.resize  (nres_0); u23_0.resize (nres_0); u24_0.resize (nres_0); u25_0.resize (nres_0);
  u26_0.resize  (nres_0); u27_0.resize  (nres_0); u28_0.resize (nres_0); u29_0.resize (nres_0); u30_0.resize (nres_0);
  u31_0.resize  (nres_0); u32_0.resize  (nres_0); u33_0.resize (nres_0); u34_0.resize (nres_0); u35_0.resize (nres_0);
  u36_0.resize  (nres_0); u37_0.resize  (nres_0); u38_0.resize (nres_0); u39_0.resize (nres_0); u40_0.resize (nres_0);
  
  if (nres_0 == nres_1)
    {
      for (int i = 0; i < nres_0; i++)
        { 
	  mres_0(i) = mres_1(i);
	  ntor_0(i) = ntor_1(i);
	}
    }
  else if (nres_0 == nres_2)
    {
      for (int i = 0; i < nres_0; i++)
	{
	  mres_0(i) = mres_2(i);
	  ntor_0(i) = ntor_2(i);
	}
    }
  else if (nres_0 == nres_3)
    {
      for (int i = 0; i < nres_0; i++)
	{
	  mres_0(i) = mres_3(i);
	  ntor_0(i) = ntor_3(i);
	}
    }
  else
    {
      for (int i = 0; i < nres_0; i++)
	{
	  mres_0(i) = mres_4(i);
	  ntor_0(i) = ntor_4(i);
	}
    }

  double weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4);
  double weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4);
  double weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4);
  double weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3);

  taua_0 = weight1 * taua_1 + weight2 * taua_2 + weight3 * taua_3 + weight4 * taua_4;
  P0_0   = weight1 * P0_1   + weight2 * P0_2   + weight3 * P0_3   + weight4 * P0_4;

  for (int j = 0; j < nres_0; j++)
    {
      u1_0 (j) = weight1 * u1_1 (j) + weight2 * u1_2 (j) + weight3 * u1_3 (j) + weight4 * u1_4 (j);
      u2_0 (j) = weight1 * u2_1 (j) + weight2 * u2_2 (j) + weight3 * u2_3 (j) + weight4 * u2_4 (j);
      u3_0 (j) = weight1 * u3_1 (j) + weight2 * u3_2 (j) + weight3 * u3_3 (j) + weight4 * u3_4 (j);
      u4_0 (j) = weight1 * u4_1 (j) + weight2 * u4_2 (j) + weight3 * u4_3 (j) + weight4 * u4_4 (j);
      u5_0 (j) = weight1 * u5_1 (j) + weight2 * u5_2 (j) + weight3 * u5_3 (j) + weight4 * u5_4 (j);
      u6_0 (j) = weight1 * u6_1 (j) + weight2 * u6_2 (j) + weight3 * u6_3 (j) + weight4 * u6_4 (j);
      u7_0 (j) = weight1 * u7_1 (j) + weight2 * u7_2 (j) + weight3 * u7_3 (j) + weight4 * u7_4 (j);
      u8_0 (j) = weight1 * u8_1 (j) + weight2 * u8_2 (j) + weight3 * u8_3 (j) + weight4 * u8_4 (j);
      u9_0 (j) = weight1 * u9_1 (j) + weight2 * u9_2 (j) + weight3 * u9_3 (j) + weight4 * u9_4 (j);
      u10_0(j) = weight1 * u10_1(j) + weight2 * u10_2(j) + weight3 * u10_3(j) + weight4 * u10_4(j);
      u11_0(j) = weight1 * u11_1(j) + weight2 * u11_2(j) + weight3 * u11_3(j) + weight4 * u11_4(j);
      u12_0(j) = weight1 * u12_1(j) + weight2 * u12_2(j) + weight3 * u12_3(j) + weight4 * u12_4(j);
      u13_0(j) = weight1 * u13_1(j) + weight2 * u13_2(j) + weight3 * u13_3(j) + weight4 * u13_4(j);
      u14_0(j) = weight1 * u14_1(j) + weight2 * u14_2(j) + weight3 * u14_3(j) + weight4 * u14_4(j);
      u15_0(j) = weight1 * u15_1(j) + weight2 * u15_2(j) + weight3 * u15_3(j) + weight4 * u15_4(j);
      u16_0(j) = weight1 * u16_1(j) + weight2 * u16_2(j) + weight3 * u16_3(j) + weight4 * u16_4(j);
      u17_0(j) = weight1 * u17_1(j) + weight2 * u17_2(j) + weight3 * u17_3(j) + weight4 * u17_4(j);
      u18_0(j) = weight1 * u18_1(j) + weight2 * u18_2(j) + weight3 * u18_3(j) + weight4 * u18_4(j);
      u19_0(j) = weight1 * u19_1(j) + weight2 * u19_2(j) + weight3 * u19_3(j) + weight4 * u19_4(j);
      u20_0(j) = weight1 * u20_1(j) + weight2 * u20_2(j) + weight3 * u20_3(j) + weight4 * u20_4(j);
      u21_0(j) = weight1 * u21_1(j) + weight2 * u21_2(j) + weight3 * u21_3(j) + weight4 * u21_4(j);
      u22_0(j) = weight1 * u22_1(j) + weight2 * u22_2(j) + weight3 * u22_3(j) + weight4 * u22_4(j);
      u23_0(j) = weight1 * u23_1(j) + weight2 * u23_2(j) + weight3 * u23_3(j) + weight4 * u23_4(j);
      u24_0(j) = weight1 * u24_1(j) + weight2 * u24_2(j) + weight3 * u24_3(j) + weight4 * u24_4(j);
      u25_0(j) = weight1 * u25_1(j) + weight2 * u25_2(j) + weight3 * u25_3(j) + weight4 * u25_4(j);
      u26_0(j) = weight1 * u26_1(j) + weight2 * u26_2(j) + weight3 * u26_3(j) + weight4 * u26_4(j);
      u27_0(j) = weight1 * u27_1(j) + weight2 * u27_2(j) + weight3 * u27_3(j) + weight4 * u27_4(j);
      u28_0(j) = weight1 * u28_1(j) + weight2 * u28_2(j) + weight3 * u28_3(j) + weight4 * u28_4(j);
      u29_0(j) = weight1 * u29_1(j) + weight2 * u29_2(j) + weight3 * u29_3(j) + weight4 * u29_4(j);
      u30_0(j) = weight1 * u30_1(j) + weight2 * u30_2(j) + weight3 * u30_3(j) + weight4 * u30_4(j);
      u31_0(j) = weight1 * u31_1(j) + weight2 * u31_2(j) + weight3 * u31_3(j) + weight4 * u31_4(j);
      u32_0(j) = weight1 * u32_1(j) + weight2 * u32_2(j) + weight3 * u32_3(j) + weight4 * u32_4(j);
      u33_0(j) = weight1 * u33_1(j) + weight2 * u33_2(j) + weight3 * u33_3(j) + weight4 * u33_4(j);
      u34_0(j) = weight1 * u34_1(j) + weight2 * u34_2(j) + weight3 * u34_3(j) + weight4 * u34_4(j);
      u35_0(j) = weight1 * u35_1(j) + weight2 * u35_2(j) + weight3 * u35_3(j) + weight4 * u35_4(j);
      u36_0(j) = weight1 * u36_1(j) + weight2 * u36_2(j) + weight3 * u36_3(j) + weight4 * u36_4(j);
      u37_0(j) = weight1 * u37_1(j) + weight2 * u37_2(j) + weight3 * u37_3(j) + weight4 * u37_4(j);
      u38_0(j) = weight1 * u38_1(j) + weight2 * u38_2(j) + weight3 * u38_3(j) + weight4 * u38_4(j);
      u39_0(j) = weight1 * u39_1(j) + weight2 * u39_2(j) + weight3 * u39_3(j) + weight4 * u39_4(j);
      u40_0(j) = weight1 * u40_1(j) + weight2 * u40_2(j) + weight3 * u40_3(j) + weight4 * u40_4(j);
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e %16.9e\n", nres_0, taua_0, P0_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j),
	     u1_0(j),  u2_0(j),  u3_0(j),  u4_0(j),  u5_0(j),
	     u6_0(j),  u7_0(j),  u8_0(j),  u9_0(j),  u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j),
	     u16_0(j), u17_0(j), u18_0(j), u19_0(j), u20_0(j),
	     u21_0(j), u22_0(j), u23_0(j), u24_0(j), u25_0(j),
	     u26_0(j), u27_0(j), u28_0(j), u29_0(j), u30_0(j),
	     u31_0(j), u32_0(j), u33_0(j), u34_0(j), u35_0(j),
	     u36_0(j), u37_0(j), u38_0(j), u39_0(j), u40_0(j));
   
  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);
  printf ("%s %11.4e\n", nFile2, weight2);
  printf ("%s %11.4e\n", nFile3, weight3);
  printf ("%s %11.4e\n", nFile4, weight4);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  fprintf (monitor, "%s %11.4e\n", nFile1, weight1);
  fprintf (monitor, "%s %11.4e\n", nFile2, weight2);
  fprintf (monitor, "%s %11.4e\n", nFile3, weight3);
  fprintf (monitor, "%s %11.4e\n", nFile4, weight4);
  fclose (monitor);
}
  
