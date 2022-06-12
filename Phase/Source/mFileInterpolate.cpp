// mFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Phase:: mFileInterp               (vector<string> mFileName,   vector<double> mFileTime,   int mFileNumber, double time)
// void Phase:: mFileInterpolateLinear    (char* mFile1, double time1, char* mFile,  double time)
// void Phase:: mFileInterpolateQuadratic (char* mFile1, double time1, char* mFile2, double time2, char* mFile,     double time)
// void Phase:: mFileInterpolateCubic     (char* mFile1, double time1, char* mFile2, double time2, char* mFile3,    double time3, char* mFile, double time)
// void Phase:: mFileInterpolateQuartic   (char* mFile1, double time1, char* mFile2, double time2, char* mFile3,    double time3,
//				           char* mFile4, double time4, char* mFile,  double time)

#include "Phase.h"

// ###############################
// Functions to interpolate mFiles
// ###############################
void Phase::mFileInterp (vector<string> mFileName, vector<double> mFileTime, int mFileNumber, double time)
{
  if (mFileNumber < 1)
    {
      printf ("PHASE::mFileInterp: Error - mFileNumber must be greater than zero\n");
      exit (1);
    }
  else if (mFileNumber == 1)
    {
      char* mFile = "Inputs/mFile";
      char* file1 = (char*) mFileName[0].c_str();
  
      mFileInterpolateLinear (file1, mFileTime[0], mFile, time);
    }
  else if (mFileNumber == 2)
    {
      char* mFile = "Inputs/mFile";
      char* file1 = (char*) mFileName[0].c_str();
      char* file2 = (char*) mFileName[1].c_str();

      mFileInterpolateQuadratic (file1, mFileTime[0], file2, mFileTime[1], mFile, time);
    }
  else if (RATS)
    {
      int index;

      if (time < mFileTime[0])
	index = 0;
      else if (time >= mFileTime[mFileNumber-1])
	index = mFileNumber - 2;
      else
	{
	  for (int i = 0; i < mFileNumber - 1; i++)
	    if (time >= mFileTime[i] && time < mFileTime[i+1])
	      index = i;
	}
      
      char* mFile = "Inputs/mFile";
      char* file1 = (char*) mFileName[index  ].c_str();
      char* file2 = (char*) mFileName[index+1].c_str();

      mFileInterpolateQuadratic (file1, mFileTime[index], file2, mFileTime[index+1], mFile, time);
    }
  else if (mFileNumber == 3)
    {
      char* mFile = "Inputs/mFile";
      char* file1 = (char*) mFileName[0].c_str();
      char* file2 = (char*) mFileName[1].c_str();
      char* file3 = (char*) mFileName[2].c_str();

      mFileInterpolateCubic (file1, mFileTime[0], file2, mFileTime[1], file3, mFileTime[2], mFile, time);
    }
  else if (mFileNumber == 4)
    {
      char* mFile = "Inputs/mFile";
      char* file1 = (char*) mFileName[0].c_str();
      char* file2 = (char*) mFileName[1].c_str();
      char* file3 = (char*) mFileName[2].c_str();
      char* file4 = (char*) mFileName[3].c_str();

      mFileInterpolateQuartic (file1, mFileTime[0], file2, mFileTime[1], file3, mFileTime[2],
			       file4, mFileTime[3], mFile, time);
    }
  else
    {
      int index, cntrl;
      
      if (time < mFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= mFileTime[mFileNumber-1])
	{
	  index = mFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < mFileNumber - 1; i++)
	    if (time >= mFileTime[i] && time < mFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == mFileNumber - 2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* mFile = "Inputs/mFile";
	  char* file1 = (char*) mFileName[index-1].c_str();
	  char* file2 = (char*) mFileName[index  ].c_str();
	  char* file3 = (char*) mFileName[index+1].c_str();
	  char* file4 = (char*) mFileName[index+2].c_str();
	  
	  mFileInterpolateQuartic (file1, mFileTime[index-1], file2, mFileTime[index], file3, mFileTime[index+1],
				   file4, mFileTime[index+2], mFile, time);
	  
	}
      else if (cntrl == 2)
	{
	  char* mFile = "Inputs/mFile";
	  char* file1 = (char*) mFileName[index  ].c_str();
	  char* file2 = (char*) mFileName[index+1].c_str();
	  char* file3 = (char*) mFileName[index+2].c_str();

	  mFileInterpolateCubic (file1, mFileTime[index], file2, mFileTime[index+1], file3, mFileTime[index+2], mFile, time);
	}
      else if (cntrl == 3)
	{
	  char* mFile = "Inputs/mFile";
	  char* file1 = (char*) mFileName[index-1].c_str();
	  char* file2 = (char*) mFileName[index  ].c_str();
	  char* file3 = (char*) mFileName[index+1].c_str();

	  mFileInterpolateCubic (file1, mFileTime[index-1], file2, mFileTime[index], file3, mFileTime[index+1], mFile, time);
	}
    }
}

void Phase::mFileInterpolateLinear (char* mFile1, double time1, char* mFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH],  line3[MAXULFILELINELENGTH], line4[MAXULFILELINELENGTH], line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineaa[MAXULFILELINELENGTH], line7[MAXULFILELINELENGTH], line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first mFile
  // ................
  FILE* file = OpenFiler (mFile1);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linea, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (lineaa, linea);

  char* token;
  token = strtok (linea, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_1 = atoi (token);

  double* v01_1 = new double[nres_1];
  double* v02_1 = new double[nres_1];
  double* v03_1 = new double[nres_1];
  double* v04_1 = new double[nres_1];
  double* v05_1 = new double[nres_1];
  double* v06_1 = new double[nres_1];
  double* v07_1 = new double[nres_1];
  double* v08_1 = new double[nres_1];
  double* v09_1 = new double[nres_1];
  double* v10_1 = new double[nres_1];
  double* v11_1 = new double[nres_1];
  double* v12_1 = new double[nres_1];
  double* v13_1 = new double[nres_1];

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i], &v13_1[i]) != 13)
      {
	printf ("PHASE:mFileInterpolateLinear: Error reading mFile_1\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate mFile data
  // ......................
  int nres = nres_1;

  double* v01 = new double[nres];
  double* v02 = new double[nres];
  double* v03 = new double[nres];
  double* v04 = new double[nres];
  double* v05 = new double[nres];
  double* v06 = new double[nres];
  double* v07 = new double[nres];
  double* v08 = new double[nres];
  double* v09 = new double[nres];
  double* v10 = new double[nres];
  double* v11 = new double[nres];
  double* v12 = new double[nres];
  double* v13 = new double[nres];

  double weight1 = 1.;

  for (int i = 0; i < nres; i++)
    {
      v01[i] = v01_1[i]; 
      v02[i] = v02_1[i]; 
      v03[i] = v03_1[i]; 
      v04[i] = v04_1[i]; 
      v05[i] = v05_1[i]; 
      v06[i] = v06_1[i]; 
      v07[i] = v07_1[i]; 
      v08[i] = v08_1[i]; 
      v09[i] = v09_1[i]; 
      v10[i] = v10_1[i]; 
      v11[i] = v11_1[i]; 
      v12[i] = v12_1[i];
      v13[i] = v13_1[i]; 
    }

  // ........................
  // Write interpolated mFile
  // ........................
  file = OpenFilew (mFile);

  fprintf (file, "%s", line1);
  fprintf (file, "%s", line2);
  fprintf (file, "%s", line3);
  fprintf (file, "%s", line4);
  fprintf (file, "%s", line5);
  fprintf (file, "%s", lineaa);
  fprintf (file, "%s", line7);
  fprintf (file, "%s", line8);

  for (int i = 0; i < nres; i++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i], v13[i]);
  
  fclose (file);

  printf ("mFile Interpolation:\n");
  printf ("%s %11.4e\n", mFile1, weight1);
 
  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1; delete[] v13_1;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12;   delete[] v13;
}

void Phase::mFileInterpolateQuadratic (char* mFile1, double time1, char* mFile2, double time2, char* mFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH], line3[MAXULFILELINELENGTH],  line4[MAXULFILELINELENGTH],  line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineb[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH], linebb[MAXULFILELINELENGTH],
       line7[MAXULFILELINELENGTH],  line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first mFile
  // ................
  FILE* file = OpenFiler (mFile1);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linea, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (lineaa, linea);

  char* token;
  token = strtok (linea, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_1 = atoi (token);

  double* v01_1 = new double[nres_1];
  double* v02_1 = new double[nres_1];
  double* v03_1 = new double[nres_1];
  double* v04_1 = new double[nres_1];
  double* v05_1 = new double[nres_1];
  double* v06_1 = new double[nres_1];
  double* v07_1 = new double[nres_1];
  double* v08_1 = new double[nres_1];
  double* v09_1 = new double[nres_1];
  double* v10_1 = new double[nres_1];
  double* v11_1 = new double[nres_1];
  double* v12_1 = new double[nres_1];
  double* v13_1 = new double[nres_1];

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i], &v13_1[i]) != 13)
      {
	printf ("PHASE:mFileInterpolateQuadratic: Error reading mFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second mFile
  // .................
  file = OpenFiler (mFile2);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (lineb, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (linebb, lineb);

  token = strtok (lineb, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_2 = atoi (token);

  double* v01_2 = new double[nres_2];
  double* v02_2 = new double[nres_2];
  double* v03_2 = new double[nres_2];
  double* v04_2 = new double[nres_2];
  double* v05_2 = new double[nres_2];
  double* v06_2 = new double[nres_2];
  double* v07_2 = new double[nres_2];
  double* v08_2 = new double[nres_2];
  double* v09_2 = new double[nres_2];
  double* v10_2 = new double[nres_2];
  double* v11_2 = new double[nres_2];
  double* v12_2 = new double[nres_2];
  double* v13_2 = new double[nres_2];
  
  for (int i = 0; i < nres_2; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_2[i], &v02_2[i], &v03_2[i], &v04_2[i], &v05_2[i], &v06_2[i], &v07_2[i], &v08_2[i], &v09_2[i], &v10_2[i], &v11_2[i], &v12_2[i], &v13_2[i]) != 13)
      {
	printf ("PHASE:mFileInterpolateQuadratic: Error reading mFile_2\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate mFile data
  // ......................
  int nres;
  if (nres_2 < nres_1)
    nres = nres_2;
  else
    nres = nres_1;
  
  int nres_0;
  if (nres_1 > nres_2)
    nres_0 = nres_1;
  else
    nres_0 = nres_2;
  
  double* v01 = new double[nres_0];
  double* v02 = new double[nres_0];
  double* v03 = new double[nres_0];
  double* v04 = new double[nres_0];
  double* v05 = new double[nres_0];
  double* v06 = new double[nres_0];
  double* v07 = new double[nres_0];
  double* v08 = new double[nres_0];
  double* v09 = new double[nres_0];
  double* v10 = new double[nres_0];
  double* v11 = new double[nres_0];
  double* v12 = new double[nres_0];
  double* v13 = new double[nres_0];

  double weight1 = (time - time2) /(time1 - time2);
  double weight2 = (time - time1) /(time2 - time1);

  for (int i = 0; i < nres; i++)
    {
      v01[i] = weight1 * v01_1[i] + weight2 * v01_2[i];
      v02[i] = weight1 * v02_1[i] + weight2 * v02_2[i];
      v03[i] = weight1 * v03_1[i] + weight2 * v03_2[i];
      v04[i] = weight1 * v04_1[i] + weight2 * v04_2[i];
      v05[i] = weight1 * v05_1[i] + weight2 * v05_2[i];
      v06[i] = weight1 * v06_1[i] + weight2 * v06_2[i];
      v07[i] = weight1 * v07_1[i] + weight2 * v07_2[i];
      v08[i] = weight1 * v08_1[i] + weight2 * v08_2[i];
      v09[i] = weight1 * v09_1[i] + weight2 * v09_2[i];
      v10[i] = weight1 * v10_1[i] + weight2 * v10_2[i];
      v11[i] = weight1 * v11_1[i] + weight2 * v11_2[i];
      v12[i] = weight1 * v12_1[i] + weight2 * v12_2[i];
      v13[i] = weight1 * v13_1[i] + weight2 * v13_2[i];
    }

  if (nres_1 > nres)
    {
      for (int i = nres; i < nres_1; i++)
	{
	  v01[i] = weight1 * v01_1[i];
	  v02[i] = weight1 * v02_1[i];
	  v03[i] = weight1 * v03_1[i];
	  v04[i] = weight1 * v04_1[i];
	  v05[i] = weight1 * v05_1[i];
	  v06[i] = weight1 * v06_1[i];
	  v07[i] = weight1 * v07_1[i];
	  v08[i] = weight1 * v08_1[i];
	  v09[i] = weight1 * v09_1[i];
	  v10[i] = weight1 * v10_1[i];
	  v11[i] = weight1 * v11_1[i];
	  v12[i] = weight1 * v12_1[i];
	  v13[i] = weight1 * v13_1[i];
	}
    }

   if (nres_2 > nres)
     {
       for (int i = nres; i < nres_2; i++)
	 {
	   v01[i] = weight2 * v01_2[i];
	   v02[i] = weight2 * v02_2[i];
	   v03[i] = weight2 * v03_2[i];
	   v04[i] = weight2 * v04_2[i];
	   v05[i] = weight2 * v05_2[i];
	   v06[i] = weight2 * v06_2[i];
	   v07[i] = weight2 * v07_2[i];
	   v08[i] = weight2 * v08_2[i];
	   v09[i] = weight2 * v09_2[i];
	   v10[i] = weight2 * v10_2[i];
	   v11[i] = weight2 * v11_2[i];
	   v12[i] = weight2 * v12_2[i];
	   v13[i] = weight2 * v13_2[i];
	 }
     }
   
  // ........................
  // Write interpolated mFile
  // ........................
  file = OpenFilew (mFile);

  fprintf (file, "%s", line1);
  fprintf (file, "%s", line2);
  fprintf (file, "%s", line3);
  fprintf (file, "%s", line4);
  fprintf (file, "%s", line5);
  if (nres_0 == nres_1)
    fprintf (file, "%s", lineaa);
  else
    fprintf (file, "%s", linebb);
  fprintf (file, "%s", line7);
  fprintf (file, "%s", line8);

  for (int i = 0; i < nres_0; i++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i], v13[i]);
  
  fclose (file);

  printf ("mFile Interpolation:\n");
  printf ("%s %11.4e\n", mFile1, weight1);
  printf ("%s %11.4e\n", mFile2, weight2);
 
  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1; delete[] v13_1;
  delete[] v01_2; delete[] v02_2; delete[] v03_2; delete[] v04_2;
  delete[] v05_2; delete[] v06_2; delete[] v07_2; delete[] v08_2;
  delete[] v09_2; delete[] v10_2; delete[] v11_2; delete[] v12_2; delete[] v13_2;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12;   delete[] v13;
}

void Phase::mFileInterpolateCubic (char* mFile1, double time1, char* mFile2, double time2, char* mFile3, double time3, char* mFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH], line3[MAXULFILELINELENGTH], line4[MAXULFILELINELENGTH],  line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineb[MAXULFILELINELENGTH], linec[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH], linebb[MAXULFILELINELENGTH],
       linecc[MAXULFILELINELENGTH], line7[MAXULFILELINELENGTH], line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first mFile
  // ................
  FILE* file = OpenFiler (mFile1);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linea, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (lineaa, linea);

  char* token;
  token = strtok (linea, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_1 = atoi (token);

  double* v01_1 = new double[nres_1];
  double* v02_1 = new double[nres_1];
  double* v03_1 = new double[nres_1];
  double* v04_1 = new double[nres_1];
  double* v05_1 = new double[nres_1];
  double* v06_1 = new double[nres_1];
  double* v07_1 = new double[nres_1];
  double* v08_1 = new double[nres_1];
  double* v09_1 = new double[nres_1];
  double* v10_1 = new double[nres_1];
  double* v11_1 = new double[nres_1];
  double* v12_1 = new double[nres_1];
  double* v13_1 = new double[nres_1];

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i], &v13_1[i]) != 13)
      {
	printf ("PHASE:mFileInterpolateCubic: Error reading mFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second mFile
  // .................
  file = OpenFiler (mFile2);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (lineb, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (linebb, lineb);

  token = strtok (lineb, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_2 = atoi (token);

  double* v01_2 = new double[nres_2];
  double* v02_2 = new double[nres_2];
  double* v03_2 = new double[nres_2];
  double* v04_2 = new double[nres_2];
  double* v05_2 = new double[nres_2];
  double* v06_2 = new double[nres_2];
  double* v07_2 = new double[nres_2];
  double* v08_2 = new double[nres_2];
  double* v09_2 = new double[nres_2];
  double* v10_2 = new double[nres_2];
  double* v11_2 = new double[nres_2];
  double* v12_2 = new double[nres_2];
  double* v13_2 = new double[nres_2];
  
  for (int i = 0; i < nres_2; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_2[i], &v02_2[i], &v03_2[i], &v04_2[i], &v05_2[i], &v06_2[i], &v07_2[i], &v08_2[i], &v09_2[i], &v10_2[i], &v11_2[i], &v12_2[i], &v13_2[i]) != 13)
      {
	printf ("PHASE:mFileInterpolateCubic: Error reading mFile_2\n");
	exit (1);
      }

  fclose (file);

  // ................
  // Read third mFile
  // ................
  file = OpenFiler (mFile3);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linec, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (linecc, linec);

  token = strtok (linec, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_3 = atoi (token);

  double* v01_3 = new double[nres_3];
  double* v02_3 = new double[nres_3];
  double* v03_3 = new double[nres_3];
  double* v04_3 = new double[nres_3];
  double* v05_3 = new double[nres_3];
  double* v06_3 = new double[nres_3];
  double* v07_3 = new double[nres_3];
  double* v08_3 = new double[nres_3];
  double* v09_3 = new double[nres_3];
  double* v10_3 = new double[nres_3];
  double* v11_3 = new double[nres_3];
  double* v12_3 = new double[nres_3];
  double* v13_3 = new double[nres_3];
  
  for (int i = 0; i < nres_3; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_3[i], &v02_3[i], &v03_3[i], &v04_3[i], &v05_3[i], &v06_3[i], &v07_3[i], &v08_3[i], &v09_3[i], &v10_3[i], &v11_3[i], &v12_3[i], &v13_3) != 13)
      {
	printf ("PHASE:mFileInterpolateCubic: Error reading mFile_3\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate mFile data
  // ......................
  if (nres_1 == nres_2 && nres_2 != nres_3)
    {
      mFileInterpolateQuadratic (mFile1, time1, mFile2, time2, mFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2)
    {
      mFileInterpolateQuadratic (mFile2, time2, mFile3, time3, mFile, time);
      return;
    }
  
  int nres;
  if (nres_2 < nres_1)
    nres = nres_2;
  else
    nres = nres_1;
  if (nres_3 < nres)
    nres = nres_3;

  double* v01 = new double[nres];
  double* v02 = new double[nres];
  double* v03 = new double[nres];
  double* v04 = new double[nres];
  double* v05 = new double[nres];
  double* v06 = new double[nres];
  double* v07 = new double[nres];
  double* v08 = new double[nres];
  double* v09 = new double[nres];
  double* v10 = new double[nres];
  double* v11 = new double[nres];
  double* v12 = new double[nres];
  double* v13 = new double[nres];

  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);

  for (int i = 0; i < nres; i++)
    {
      v01[i] = weight1 * v01_1[i] + weight2 * v01_2[i] + weight3 * v01_3[i];
      v02[i] = weight1 * v02_1[i] + weight2 * v02_2[i] + weight3 * v02_3[i];
      v03[i] = weight1 * v03_1[i] + weight2 * v03_2[i] + weight3 * v03_3[i];
      v04[i] = weight1 * v04_1[i] + weight2 * v04_2[i] + weight3 * v04_3[i];
      v05[i] = weight1 * v05_1[i] + weight2 * v05_2[i] + weight3 * v05_3[i];
      v06[i] = weight1 * v06_1[i] + weight2 * v06_2[i] + weight3 * v06_3[i];
      v07[i] = weight1 * v07_1[i] + weight2 * v07_2[i] + weight3 * v07_3[i];
      v08[i] = weight1 * v08_1[i] + weight2 * v08_2[i] + weight3 * v08_3[i];
      v09[i] = weight1 * v09_1[i] + weight2 * v09_2[i] + weight3 * v09_3[i];
      v10[i] = weight1 * v10_1[i] + weight2 * v10_2[i] + weight3 * v10_3[i];
      v11[i] = weight1 * v11_1[i] + weight2 * v11_2[i] + weight3 * v11_3[i];
      v12[i] = weight1 * v12_1[i] + weight2 * v12_2[i] + weight3 * v12_3[i];
      v13[i] = weight1 * v13_1[i] + weight2 * v13_2[i] + weight3 * v13_3[i];
    }

  // ........................
  // Write interpolated mFile
  // ........................
  file = OpenFilew (mFile);

  fprintf (file, "%s", line1);
  fprintf (file, "%s", line2);
  fprintf (file, "%s", line3);
  fprintf (file, "%s", line4);
  fprintf (file, "%s", line5);
  if (nres == nres_1)
    fprintf (file, "%s", lineaa);
  else if (nres == nres_2)
    fprintf (file, "%s", linebb);
  else
    fprintf (file, "%s", linecc);
  fprintf (file, "%s", line7);
  fprintf (file, "%s", line8);

  for (int i = 0; i < nres; i++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i], v13[i]);
  
  fclose (file);

  printf ("mFile Interpolation:\n");
  printf ("%s %11.4e\n", mFile1, weight1);
  printf ("%s %11.4e\n", mFile2, weight2);
  printf ("%s %11.4e\n", mFile3, weight3);

  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1; delete[] v13_1;
  delete[] v01_2; delete[] v02_2; delete[] v03_2; delete[] v04_2;
  delete[] v05_2; delete[] v06_2; delete[] v07_2; delete[] v08_2;
  delete[] v09_2; delete[] v10_2; delete[] v11_2; delete[] v12_2; delete[] v13_2;
  delete[] v01_3; delete[] v02_3; delete[] v03_3; delete[] v04_3;
  delete[] v05_3; delete[] v06_3; delete[] v07_3; delete[] v08_3;
  delete[] v09_3; delete[] v10_3; delete[] v11_3; delete[] v12_3; delete[] v13_3;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12;   delete[] v13; 
}

void Phase::mFileInterpolateQuartic (char* mFile1, double time1, char* mFile2, double time2, char* mFile3, double time3,
				     char* mFile4, double time4, char* mFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH],  line3[MAXULFILELINELENGTH],  line4[MAXULFILELINELENGTH], line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineb[MAXULFILELINELENGTH],  linec[MAXULFILELINELENGTH],  lined[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH],
       linebb[MAXULFILELINELENGTH], linecc[MAXULFILELINELENGTH], linedd[MAXULFILELINELENGTH], line7[MAXULFILELINELENGTH], line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first mFile
  // ................
  FILE* file = OpenFiler (mFile1);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linea, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (lineaa, linea);

  char* token;
  token = strtok (linea, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_1 = atoi (token);

  double* v01_1 = new double[nres_1];
  double* v02_1 = new double[nres_1];
  double* v03_1 = new double[nres_1];
  double* v04_1 = new double[nres_1];
  double* v05_1 = new double[nres_1];
  double* v06_1 = new double[nres_1];
  double* v07_1 = new double[nres_1];
  double* v08_1 = new double[nres_1];
  double* v09_1 = new double[nres_1];
  double* v10_1 = new double[nres_1];
  double* v11_1 = new double[nres_1];
  double* v12_1 = new double[nres_1];
  double* v13_1 = new double[nres_1];

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i], &v13_1[i]) != 13)
      {
	printf ("PHASE::mFileInterpolateQuartic: Error reading mFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second mFile
  // .................
  file = OpenFiler (mFile2);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (lineb, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (linebb, lineb);

  token = strtok (lineb, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_2 = atoi (token);

  double* v01_2 = new double[nres_2];
  double* v02_2 = new double[nres_2];
  double* v03_2 = new double[nres_2];
  double* v04_2 = new double[nres_2];
  double* v05_2 = new double[nres_2];
  double* v06_2 = new double[nres_2];
  double* v07_2 = new double[nres_2];
  double* v08_2 = new double[nres_2];
  double* v09_2 = new double[nres_2];
  double* v10_2 = new double[nres_2];
  double* v11_2 = new double[nres_2];
  double* v12_2 = new double[nres_2];
  double* v13_2 = new double[nres_2];
  
  for (int i = 0; i < nres_2; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_2[i], &v02_2[i], &v03_2[i], &v04_2[i], &v05_2[i], &v06_2[i], &v07_2[i], &v08_2[i], &v09_2[i], &v10_2[i], &v11_2[i], &v12_2[i], &v13_2[i]) != 13)
      {
	printf ("PHASE::mFileInterpolateQuartic: Error reading mFile_2\n");
	exit (1);
      }

  fclose (file);

  // ................
  // Read third mFile
  // ................
  file = OpenFiler (mFile3);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linec, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (linecc, linec);

  token = strtok (linec, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_3 = atoi (token);

  double* v01_3 = new double[nres_3];
  double* v02_3 = new double[nres_3];
  double* v03_3 = new double[nres_3];
  double* v04_3 = new double[nres_3];
  double* v05_3 = new double[nres_3];
  double* v06_3 = new double[nres_3];
  double* v07_3 = new double[nres_3];
  double* v08_3 = new double[nres_3];
  double* v09_3 = new double[nres_3];
  double* v10_3 = new double[nres_3];
  double* v11_3 = new double[nres_3];
  double* v12_3 = new double[nres_3];
  double* v13_3 = new double[nres_3];
  
  for (int i = 0; i < nres_3; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_3[i], &v02_3[i], &v03_3[i], &v04_3[i], &v05_3[i], &v06_3[i], &v07_3[i], &v08_3[i], &v09_3[i], &v10_3[i], &v11_3[i], &v12_3[i], &v13_3[i]) != 13)
      {
	printf ("PHASE::mFileInterpolateQuartic: Error reading mFile_3\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read fourth mFile
  // .................
  file = OpenFiler (mFile4);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (lined, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (linedd, lined);

  token = strtok (lined, " "); 
  token = strtok (NULL,  " "); 
  token = strtok (NULL,  " ");
  int nres_4 = atoi (token);

  double* v01_4 = new double[nres_4];
  double* v02_4 = new double[nres_4];
  double* v03_4 = new double[nres_4];
  double* v04_4 = new double[nres_4];
  double* v05_4 = new double[nres_4];
  double* v06_4 = new double[nres_4];
  double* v07_4 = new double[nres_4];
  double* v08_4 = new double[nres_4];
  double* v09_4 = new double[nres_4];
  double* v10_4 = new double[nres_4];
  double* v11_4 = new double[nres_4];
  double* v12_4 = new double[nres_4];
  double* v13_4 = new double[nres_4];
  
  for (int i = 0; i < nres_4; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_4[i], &v02_4[i], &v03_4[i], &v04_4[i], &v05_4[i], &v06_4[i], &v07_4[i], &v08_4[i], &v09_4[i], &v10_4[i], &v11_4[i], &v12_4[i], &v13_4[i]) != 13)
      {
	printf ("PHASE::mFileInterpolateQuartic: Error reading mFile_4\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate mFile data
  // ......................
  if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 != nres_2)
    {
      mFileInterpolateQuadratic (mFile2, time2, mFile3, time3, mFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 == nres_2 && nres_4 != nres_2)
    {
      mFileInterpolateCubic (mFile1, time1, mFile2, time2, mFile3, time3, mFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 == nres_2)
    {
      mFileInterpolateCubic (mFile2, time2, mFile3, time3, mFile4, time4, mFile, time);
      return;
    }
   
  int nres;
  if (nres_2 < nres_1)
    nres = nres_2;
  else
    nres = nres_1;
  if (nres_3 < nres)
    nres = nres_3;
  if (nres_4 < nres)
    nres = nres_4;

  double* v01 = new double[nres];
  double* v02 = new double[nres];
  double* v03 = new double[nres];
  double* v04 = new double[nres];
  double* v05 = new double[nres];
  double* v06 = new double[nres];
  double* v07 = new double[nres];
  double* v08 = new double[nres];
  double* v09 = new double[nres];
  double* v10 = new double[nres];
  double* v11 = new double[nres];
  double* v12 = new double[nres];
  double* v13 = new double[nres];

  double weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4);
  double weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4);
  double weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4);
  double weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3);

  for (int i = 0; i < nres; i++)
    {
      v01[i] = weight1 * v01_1[i] + weight2 * v01_2[i] + weight3 * v01_3[i] + weight4 * v01_4[i];
      v02[i] = weight1 * v02_1[i] + weight2 * v02_2[i] + weight3 * v02_3[i] + weight4 * v02_4[i];
      v03[i] = weight1 * v03_1[i] + weight2 * v03_2[i] + weight3 * v03_3[i] + weight4 * v03_4[i];
      v04[i] = weight1 * v04_1[i] + weight2 * v04_2[i] + weight3 * v04_3[i] + weight4 * v04_4[i];
      v05[i] = weight1 * v05_1[i] + weight2 * v05_2[i] + weight3 * v05_3[i] + weight4 * v05_4[i];
      v06[i] = weight1 * v06_1[i] + weight2 * v06_2[i] + weight3 * v06_3[i] + weight4 * v06_4[i];
      v07[i] = weight1 * v07_1[i] + weight2 * v07_2[i] + weight3 * v07_3[i] + weight4 * v07_4[i];
      v08[i] = weight1 * v08_1[i] + weight2 * v08_2[i] + weight3 * v08_3[i] + weight4 * v08_4[i];
      v09[i] = weight1 * v09_1[i] + weight2 * v09_2[i] + weight3 * v09_3[i] + weight4 * v09_4[i];
      v10[i] = weight1 * v10_1[i] + weight2 * v10_2[i] + weight3 * v10_3[i] + weight4 * v10_4[i];
      v11[i] = weight1 * v11_1[i] + weight2 * v11_2[i] + weight3 * v11_3[i] + weight4 * v11_4[i];
      v12[i] = weight1 * v12_1[i] + weight2 * v12_2[i] + weight3 * v12_3[i] + weight4 * v12_4[i];
      v13[i] = weight1 * v13_1[i] + weight2 * v13_2[i] + weight3 * v13_3[i] + weight4 * v13_4[i];
    }

  // ........................
  // Write interpolated mFile
  // ........................
  file = OpenFilew (mFile);

  fprintf (file, "%s", line1);
  fprintf (file, "%s", line2);
  fprintf (file, "%s", line3);
  fprintf (file, "%s", line4);
  fprintf (file, "%s", line5);
  if (nres == nres_1)
    fprintf (file, "%s", lineaa);
  else if (nres == nres_2)
    fprintf (file, "%s", linebb);
  else if (nres == nres_3)
    fprintf (file, "%s", linecc);
  else
    fprintf (file, "%s", linedd);
  fprintf (file, "%s", line7);
  fprintf (file, "%s", line8);

  for (int i = 0; i < nres; i++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i], v13[i]);
  
  fclose (file);

  printf ("mFile Interpolation:\n");
  printf ("%s %11.4e\n", mFile1, weight1);
  printf ("%s %11.4e\n", mFile2, weight2);
  printf ("%s %11.4e\n", mFile3, weight3);
  printf ("%s %11.4e\n", mFile4, weight4);

  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1; delete[] v13_1;
  delete[] v01_2; delete[] v02_2; delete[] v03_2; delete[] v04_2;
  delete[] v05_2; delete[] v06_2; delete[] v07_2; delete[] v08_2;
  delete[] v09_2; delete[] v10_2; delete[] v11_2; delete[] v12_2; delete[] v13_2;
  delete[] v01_3; delete[] v02_3; delete[] v03_3; delete[] v04_3;
  delete[] v05_3; delete[] v06_3; delete[] v07_3; delete[] v08_3;
  delete[] v09_3; delete[] v10_3; delete[] v11_3; delete[] v12_3; delete[] v13_3;
  delete[] v01_4; delete[] v02_4; delete[] v03_4; delete[] v04_4;
  delete[] v05_4; delete[] v06_4; delete[] v07_4; delete[] v08_4;
  delete[] v09_4; delete[] v10_4; delete[] v11_4; delete[] v12_4; delete[] v13_4;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12;   delete[] v13;
}

