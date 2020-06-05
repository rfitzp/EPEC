// nFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Phase:: nFileInterp               (vector<string> nFileName,   vector<double> nFileTime,   int nFileNumber, double time)
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
  if (nFileNumber < 2)
    {
      printf ("PHASE::nFileInterp: Error - nFileNumber must be greater than unity\n");
      exit (1);
    }
  else if (nFileNumber == 2)
    {
      char* nFile = "../Neoclassical/nFile";
      char* file1 = (char*) nFileName[0].c_str();
      char* file2 = (char*) nFileName[1].c_str();

      nFileInterpolateQuadratic (file1, nFileTime[0], file2, nFileTime[1], nFile, time);
    }
  else if (nFileNumber == 3)
    {
      char* nFile = "../Neoclassical/nFile";
      char* file1 = (char*) nFileName[0].c_str();
      char* file2 = (char*) nFileName[1].c_str();
      char* file3 = (char*) nFileName[2].c_str();

      nFileInterpolateCubic (file1, nFileTime[0], file2, nFileTime[1], file3, nFileTime[2], nFile, time);
    }
  else if (nFileNumber == 4)
    {
      char* nFile = "../Neoclassical/nFile";
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
	  char* nFile = "../Neoclassical/nFile";
	  char* file1 = (char*) nFileName[index-1].c_str();
	  char* file2 = (char*) nFileName[index  ].c_str();
	  char* file3 = (char*) nFileName[index+1].c_str();
	  char* file4 = (char*) nFileName[index+2].c_str();
	  
	  nFileInterpolateQuartic (file1, nFileTime[index-1], file2, nFileTime[index], file3, nFileTime[index+1],
				   file4, nFileTime[index+2], nFile, time);
	}
      else if (cntrl == 2)
	{
	  char* nFile = "../Neoclassical/nFile";
	  char* file1 = (char*) nFileName[index  ].c_str();
	  char* file2 = (char*) nFileName[index+1].c_str();
	  char* file3 = (char*) nFileName[index+2].c_str();
	  
	  nFileInterpolateCubic (file1, nFileTime[index], file2, nFileTime[index+1], file3, nFileTime[index+2], nFile, time);
	}
      else if (cntrl == 3)
	{
	  char* nFile = "../Neoclassical/nFile";
	  char* file1 = (char*) nFileName[index-1].c_str();
	  char* file2 = (char*) nFileName[index  ].c_str();
	  char* file3 = (char*) nFileName[index+1].c_str();
	  
	  nFileInterpolateCubic (file1, nFileTime[index-1], file2, nFileTime[index], file3, nFileTime[index+1], nFile, time);
	}
    }
}

void Phase::nFileInterpolateQuadratic (char* nFile1, double time1, char* nFile2, double time2, char* nFile, double time)
{
  int             nres_1;
  double          taua_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;

  int             nres_2;
  double          taua_2;
  Array<int, 1>   mres_2, ntor_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2;
  Array<double,1> u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2, u18_2, u19_2;

  int             nres_0;
  double          taua_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
 
  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf", &nres_1, &taua_1) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1);
  ntor_1.resize (nres_1);
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
  u11_1.resize  (nres_1);
  u12_1.resize  (nres_1);
  u13_1.resize  (nres_1);
  u14_1.resize  (nres_1);
  u15_1.resize  (nres_1);
  u16_1.resize  (nres_1);
  u17_1.resize  (nres_1);
  u18_1.resize  (nres_1);
  u19_1.resize  (nres_1);
  
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // .................
  // Read second nFile
  // .................
  file = OpenFiler (nFile2);
  
  if (fscanf (file, "%d %lf", &nres_2, &taua_2) != 2)
    {
      printf ("PHASE::nFileInterpolateQuadratic:Error reading nFile_2 (1)\n");
      exit (1);
    }
  
  mres_2.resize (nres_2);
  ntor_2.resize (nres_2);
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
  u11_2.resize  (nres_2);
  u12_2.resize  (nres_2);
  u13_2.resize  (nres_2);
  u14_2.resize  (nres_2);
  u15_2.resize  (nres_2);
  u16_2.resize  (nres_2);
  u17_2.resize  (nres_2);
  u18_2.resize  (nres_2);
  u19_2.resize  (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &ntor_2(j),
		  &u1_2 (j), &u2_2 (j), &u3_2 (j), &u4_2 (j), &u5_2 (j),
		  &u6_2 (j), &u7_2 (j), &u8_2 (j), &u9_2 (j), &u10_2(j),
		  &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j),
		  &u16_2(j), &u17_2(j), &u18_2(j), &u19_2(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_2 (2)\n");
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

  mres_0.resize(nres_0);
  ntor_0.resize(nres_0);
  u1_0.resize  (nres_0);
  u2_0.resize  (nres_0);
  u3_0.resize  (nres_0);
  u4_0.resize  (nres_0);
  u5_0.resize  (nres_0);
  u6_0.resize  (nres_0);
  u7_0.resize  (nres_0);
  u8_0.resize  (nres_0);
  u9_0.resize  (nres_0);
  u10_0.resize (nres_0);
  u11_0.resize (nres_0);
  u12_0.resize (nres_0);
  u13_0.resize (nres_0);
  u14_0.resize (nres_0);
  u15_0.resize (nres_0);
  u16_0.resize (nres_0);
  u17_0.resize (nres_0);
  u18_0.resize (nres_0);
  u19_0.resize (nres_0);

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
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e\n", nres_0, taua_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j), u18_0(j), u19_0(j));

  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);
  printf ("%s %11.4e\n", nFile2, weight2);
}
  
void Phase::nFileInterpolateCubic (char* nFile1, double time1, char* nFile2, double time2, char* nFile3, double time3, char* nFile, double time)
{
  int             nres_1;
  double          taua_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;

  int             nres_2;
  double          taua_2;
  Array<int, 1>   mres_2, ntor_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2;
  Array<double,1> u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2, u18_2, u19_2;

  int             nres_3;
  double          taua_3;
  Array<int, 1>   mres_3, ntor_3;
  Array<double,1> u1_3, u2_3, u3_3, u4_3, u5_3, u6_3, u7_3, u8_3, u9_3, u10_3;
  Array<double,1> u11_3, u12_3, u13_3, u14_3, u15_3, u16_3, u17_3, u18_3, u19_3;

  int             nres_0;
  double          taua_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
 
  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf", &nres_1, &taua_1) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1);
  ntor_1.resize (nres_1);
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
  u11_1.resize  (nres_1);
  u12_1.resize  (nres_1);
  u13_1.resize  (nres_1);
  u14_1.resize  (nres_1);
  u15_1.resize  (nres_1);
  u16_1.resize  (nres_1);
  u17_1.resize  (nres_1);
  u18_1.resize  (nres_1);
  u19_1.resize  (nres_1);
  
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // .................
  // Read second nFile
  // .................
  file = OpenFiler (nFile2);
  
  if (fscanf (file, "%d %lf", &nres_2, &taua_2) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_2 (1)\n");
      exit (1);
    }
  
  mres_2.resize (nres_2);
  ntor_2.resize (nres_2);
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
  u11_2.resize  (nres_2);
  u12_2.resize  (nres_2);
  u13_2.resize  (nres_2);
  u14_2.resize  (nres_2);
  u15_2.resize  (nres_2);
  u16_2.resize  (nres_2);
  u17_2.resize  (nres_2);
  u18_2.resize  (nres_2);
  u19_2.resize  (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &ntor_2(j),
		  &u1_2 (j), &u2_2 (j), &u3_2 (j), &u4_2 (j), &u5_2 (j),
		  &u6_2 (j), &u7_2 (j), &u8_2 (j), &u9_2 (j), &u10_2(j),
		  &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j),
		  &u16_2(j), &u17_2(j), &u18_2(j), &u19_2(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_2 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ................
  // Read third nFile
  // ................
  file = OpenFiler (nFile3);

  if (fscanf (file, "%d %lf", &nres_3, &taua_3) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_3 (1)\n");
      exit (1);
    }

  mres_3.resize (nres_3);
  ntor_3.resize (nres_3);
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
  u11_3.resize  (nres_3);
  u12_3.resize  (nres_3);
  u13_3.resize  (nres_3);
  u14_3.resize  (nres_3);
  u15_3.resize  (nres_3);
  u16_3.resize  (nres_3);
  u17_3.resize  (nres_3);
  u18_3.resize  (nres_3);
  u19_3.resize  (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &ntor_3(j),
		  &u1_3 (j), &u2_3 (j), &u3_3 (j), &u4_3 (j), &u5_3 (j),
		  &u6_3 (j), &u7_3 (j), &u8_3 (j), &u9_3 (j), &u10_3(j),
		  &u11_3(j), &u12_3(j), &u13_3(j), &u14_3(j), &u15_3(j),
		  &u16_3(j), &u17_3(j), &u18_3(j), &u19_3(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_3 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ......................
  // Interpolate nFile data
  // ......................
  if (nres_1 == nres_2 && nres_2 != nres_3)
    nFileInterpolateQuadratic (nFile1, time1, nFile2, time2, nFile, time);
  else if (nres_2 == nres_3 && nres_1 != nres_2)
    nFileInterpolateQuadratic (nFile2, time2, nFile3, time3, nFile, time);
  
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
  if (nres_3 < nres_0)
    nres_0 = nres_3;
  
  mres_0.resize(nres_0);
  ntor_0.resize(nres_0);
  u1_0.resize  (nres_0);
  u2_0.resize  (nres_0);
  u3_0.resize  (nres_0);
  u4_0.resize  (nres_0);
  u5_0.resize  (nres_0);
  u6_0.resize  (nres_0);
  u7_0.resize  (nres_0);
  u8_0.resize  (nres_0);
  u9_0.resize  (nres_0);
  u10_0.resize (nres_0);
  u11_0.resize (nres_0);
  u12_0.resize (nres_0);
  u13_0.resize (nres_0);
  u14_0.resize (nres_0);
  u15_0.resize (nres_0);
  u16_0.resize (nres_0);
  u17_0.resize (nres_0);
  u18_0.resize (nres_0);
  u19_0.resize (nres_0);

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
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e\n", nres_0, taua_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j), u18_0(j), u19_0(j));

  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);
  printf ("%s %11.4e\n", nFile2, weight2);
  printf ("%s %11.4e\n", nFile3, weight3);
}

void Phase::nFileInterpolateQuartic (char* nFile1, double time1, char* nFile2, double time2, char* nFile3, double time3,
				     char* nFile4, double time4, char* nFile, double time)
{
  int             nres_1;
  double          taua_1;
  Array<int, 1>   mres_1, ntor_1;
  Array<double,1> u1_1, u2_1, u3_1, u4_1, u5_1, u6_1, u7_1, u8_1, u9_1, u10_1;
  Array<double,1> u11_1, u12_1, u13_1, u14_1, u15_1, u16_1, u17_1, u18_1, u19_1;

  int             nres_2;
  double          taua_2;
  Array<int, 1>   mres_2, ntor_2;
  Array<double,1> u1_2, u2_2, u3_2, u4_2, u5_2, u6_2, u7_2, u8_2, u9_2, u10_2;
  Array<double,1> u11_2, u12_2, u13_2, u14_2, u15_2, u16_2, u17_2, u18_2, u19_2;

  int             nres_3;
  double          taua_3;
  Array<int, 1>   mres_3, ntor_3;
  Array<double,1> u1_3, u2_3, u3_3, u4_3, u5_3, u6_3, u7_3, u8_3, u9_3, u10_3;
  Array<double,1> u11_3, u12_3, u13_3, u14_3, u15_3, u16_3, u17_3, u18_3, u19_3;

  int             nres_4;
  double          taua_4;
  Array<int, 1>   mres_4, ntor_4;
  Array<double,1> u1_4, u2_4, u3_4, u4_4, u5_4, u6_4, u7_4, u8_4, u9_4, u10_4;
  Array<double,1> u11_4, u12_4, u13_4, u14_4, u15_4, u16_4, u17_4, u18_4, u19_4;
  
  int             nres_0;
  double          taua_0;
  Array<int, 1>   mres_0, ntor_0;
  Array<double,1> u1_0, u2_0, u3_0, u4_0, u5_0, u6_0, u7_0, u8_0, u9_0, u10_0;
  Array<double,1> u11_0, u12_0, u13_0, u14_0, u15_0, u16_0, u17_0, u18_0, u19_0;
 
  // ................
  // Read first nFile
  // ................
  FILE* file = OpenFiler (nFile1);
 
  if (fscanf (file, "%d %lf", &nres_1, &taua_1) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_1 (1)\n");
      exit (1);
    }
  
  mres_1.resize (nres_1);
  ntor_1.resize (nres_1);
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
  u11_1.resize  (nres_1);
  u12_1.resize  (nres_1);
  u13_1.resize  (nres_1);
  u14_1.resize  (nres_1);
  u15_1.resize  (nres_1);
  u16_1.resize  (nres_1);
  u17_1.resize  (nres_1);
  u18_1.resize  (nres_1);
  u19_1.resize  (nres_1);
  
  for (int j = 0; j < nres_1; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_1(j), &ntor_1(j),
		  &u1_1 (j), &u2_1 (j), &u3_1 (j), &u4_1 (j), &u5_1 (j),
		  &u6_1 (j), &u7_1 (j), &u8_1 (j), &u9_1 (j), &u10_1(j),
		  &u11_1(j), &u12_1(j), &u13_1(j), &u14_1(j), &u15_1(j),
		  &u16_1(j), &u17_1(j), &u18_1(j), &u19_1(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_1 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // .................
  // Read second nFile
  // .................
  file = OpenFiler (nFile2);
  
  if (fscanf (file, "%d %lf", &nres_2, &taua_2) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_2 (1)\n");
      exit (1);
    }
  
  mres_2.resize (nres_2);
  ntor_2.resize (nres_2);
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
  u11_2.resize  (nres_2);
  u12_2.resize  (nres_2);
  u13_2.resize  (nres_2);
  u14_2.resize  (nres_2);
  u15_2.resize  (nres_2);
  u16_2.resize  (nres_2);
  u17_2.resize  (nres_2);
  u18_2.resize  (nres_2);
  u19_2.resize  (nres_2);
  
  for (int j = 0; j < nres_2; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_2(j), &ntor_2(j),
		  &u1_2 (j), &u2_2 (j), &u3_2 (j), &u4_2 (j), &u5_2 (j),
		  &u6_2 (j), &u7_2 (j), &u8_2 (j), &u9_2 (j), &u10_2(j),
		  &u11_2(j), &u12_2(j), &u13_2(j), &u14_2(j), &u15_2(j),
		  &u16_2(j), &u17_2(j), &u18_2(j), &u19_2(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_2 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ................
  // Read third nFile
  // ................
  file = OpenFiler (nFile3);

  if (fscanf (file, "%d %lf", &nres_3, &taua_3) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_3 (1)\n");
      exit (1);
    }

  mres_3.resize (nres_3);
  ntor_3.resize (nres_3);
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
  u11_3.resize  (nres_3);
  u12_3.resize  (nres_3);
  u13_3.resize  (nres_3);
  u14_3.resize  (nres_3);
  u15_3.resize  (nres_3);
  u16_3.resize  (nres_3);
  u17_3.resize  (nres_3);
  u18_3.resize  (nres_3);
  u19_3.resize  (nres_3);
  
  for (int j = 0; j < nres_3; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_3(j), &ntor_3(j),
		  &u1_3 (j), &u2_3 (j), &u3_3 (j), &u4_3 (j), &u5_3 (j),
		  &u6_3 (j), &u7_3 (j), &u8_3 (j), &u9_3 (j), &u10_3(j),
		  &u11_3(j), &u12_3(j), &u13_3(j), &u14_3(j), &u15_3(j),
		  &u16_3(j), &u17_3(j), &u18_3(j), &u19_3(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_3 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);

  // .................
  // Read fourth nFile
  // .................
  file = OpenFiler (nFile4);

  if (fscanf (file, "%d %lf", &nres_4, &taua_4) != 2)
    {
      printf ("PHASE::nFileInterpolateCubic:Error reading nFile_4 (1)\n");
      exit (1);
    }

  mres_4.resize (nres_4);
  ntor_4.resize (nres_4);
  u1_4.resize   (nres_4);
  u2_4.resize   (nres_4);
  u3_4.resize   (nres_4);
  u4_4.resize   (nres_4);
  u5_4.resize   (nres_4);
  u6_4.resize   (nres_4);
  u7_4.resize   (nres_4);
  u8_4.resize   (nres_4);
  u9_4.resize   (nres_4);
  u10_4.resize  (nres_4);
  u11_4.resize  (nres_4);
  u12_4.resize  (nres_4);
  u13_4.resize  (nres_4);
  u14_4.resize  (nres_4);
  u15_4.resize  (nres_4);
  u16_4.resize  (nres_4);
  u17_4.resize  (nres_4);
  u18_4.resize  (nres_4);
  u19_4.resize  (nres_4);
  
  for (int j = 0; j < nres_4; j++)
    {
      if (fscanf (file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &mres_4(j), &ntor_4(j),
		  &u1_4 (j), &u2_4 (j), &u3_4 (j), &u4_4 (j), &u5_4 (j),
		  &u6_4 (j), &u7_4 (j), &u8_4 (j), &u9_4 (j), &u10_4(j),
		  &u11_4(j), &u12_4(j), &u13_4(j), &u14_4(j), &u15_4(j),
		  &u16_4(j), &u17_4(j), &u18_4(j), &u19_4(j)) != 21)
	{
	  printf ("NEOCLASSICAL: Error reading nFile_4 (2)\n");
	  exit (1);
	}
    }
  
  fclose (file);
  
  // ......................
  // Interpolate nFile data
  // ......................
  if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 != nres_2)
    nFileInterpolateQuadratic (nFile2, time2, nFile3, time3, nFile, time);
  else if (nres_2 == nres_3 && nres_1 == nres_2 && nres_4 != nres_2)
    nFileInterpolateCubic (nFile1, time1, nFile2, time2, nFile3, time3, nFile, time);
  else if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 == nres_2)
    nFileInterpolateCubic (nFile2, time2, nFile3, time3, nFile4, time4, nFile, time);
			   
  if (nres_2 < nres_1)
    nres_0 = nres_2;
  else
    nres_0 = nres_1;
  if (nres_3 < nres_0)
    nres_0 = nres_3;
  if (nres_4 < nres_0)
    nres_0 = nres_4;
  
  mres_0.resize(nres_0);
  ntor_0.resize(nres_0);
  u1_0.resize  (nres_0);
  u2_0.resize  (nres_0);
  u3_0.resize  (nres_0);
  u4_0.resize  (nres_0);
  u5_0.resize  (nres_0);
  u6_0.resize  (nres_0);
  u7_0.resize  (nres_0);
  u8_0.resize  (nres_0);
  u9_0.resize  (nres_0);
  u10_0.resize (nres_0);
  u11_0.resize (nres_0);
  u12_0.resize (nres_0);
  u13_0.resize (nres_0);
  u14_0.resize (nres_0);
  u15_0.resize (nres_0);
  u16_0.resize (nres_0);
  u17_0.resize (nres_0);
  u18_0.resize (nres_0);
  u19_0.resize (nres_0);

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
    }
  
  // ........................
  // Write interpolated nFile
  // ........................
  file = OpenFilew (nFile);

  fprintf (file, "%3d %16.9e\n", nres_0, taua_0);
 
  for (int j = 0; j < nres_0; j++)
    fprintf (file, "%d %d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres_0(j), ntor_0(j), u1_0(j), u2_0(j), u3_0(j), u4_0(j), u5_0(j), u6_0(j), u7_0(j), u8_0(j), u9_0(j), u10_0(j),
	     u11_0(j), u12_0(j), u13_0(j), u14_0(j), u15_0(j), u16_0(j), u17_0(j), u18_0(j), u19_0(j));

  fclose (file);
  
  printf ("nFile Interpolation:\n");
  printf ("%s %11.4e\n", nFile1, weight1);
  printf ("%s %11.4e\n", nFile2, weight2);
  printf ("%s %11.4e\n", nFile3, weight3);
  printf ("%s %11.4e\n", nFile4, weight4);
}
  
