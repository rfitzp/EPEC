// lFileInterpolate.h

// PROGRAM ORGANIZATION:
//
// void Phase:: lFileInterp               (vector<string> lFileName,   vector<double> lFileTime,   int lFileNumber, double time)
// void Phase:: lFileInterpolateQuadratic (char* lFile1, double time1, char* lFile2, double time2, char* lFile,     double time)
// void Phase:: lFileInterpolateCubic     (char* lFile1, double time1, char* lFile2, double time2, char* lFile3,    double time3, char* lFile, double time)
// void Phase:: lFileInterpolateQuartic   (char* lFile1, double time1, char* lFile2, double time2, char* lFile3,    double time3,
//				           char* lFile4, double time4, char* lFile,  double time)

#include "Phase.h"

// ###############################
// Functions to interpolate lFiles
// ###############################
void Phase::lFileInterp (vector<string> lFileName, vector<double> lFileTime, int lFileNumber, double time)
{
  if (lFileNumber < 2)
    {
      printf ("PHASE::lFileInterp: Error - lFileNumber must be greater than unity\n");
      exit (1);
    }
  else if (lFileNumber == 2)
    {
      char* lFile = "Inputs/lFile";
      char* file1 = (char*) lFileName[0].c_str();
      char* file2 = (char*) lFileName[1].c_str();

      lFileInterpolateQuadratic (file1, lFileTime[0], file2, lFileTime[1], lFile, time);
    }
  else if (lFileNumber == 3)
    {
      char* lFile = "Inputs/lFile";
      char* file1 = (char*) lFileName[0].c_str();
      char* file2 = (char*) lFileName[1].c_str();
      char* file3 = (char*) lFileName[2].c_str();

      lFileInterpolateCubic (file1, lFileTime[0], file2, lFileTime[1], file3, lFileTime[2], lFile, time);
    }
  else if (lFileNumber == 4)
    {
      char* lFile = "Inputs/lFile";
      char* file1 = (char*) lFileName[0].c_str();
      char* file2 = (char*) lFileName[1].c_str();
      char* file3 = (char*) lFileName[2].c_str();
      char* file4 = (char*) lFileName[3].c_str();

      lFileInterpolateQuartic (file1, lFileTime[0], file2, lFileTime[1], file3, lFileTime[2],
			       file4, lFileTime[3], lFile, time);
    }
  else
    {
      int index, cntrl;

      if (time < lFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= lFileTime[lFileNumber-1])
	{
	  index = lFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < lFileNumber-1; i++)
	    if (time >= lFileTime[i] && time < lFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == lFileNumber-2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* lFile = "Inputs/lFile";
	  char* file1 = (char*) lFileName[index-1].c_str();
	  char* file2 = (char*) lFileName[index  ].c_str();
	  char* file3 = (char*) lFileName[index+1].c_str();
	  char* file4 = (char*) lFileName[index+2].c_str();
	  
	  lFileInterpolateQuartic (file1, lFileTime[index-1], file2, lFileTime[index], file3, lFileTime[index+1],
				   file4, lFileTime[index+2], lFile, time);
	}
      else if (cntrl == 2)
	{
	  char* lFile = "Inputs/lFile";
	  char* file1 = (char*) lFileName[index  ].c_str();
	  char* file2 = (char*) lFileName[index+1].c_str();
	  char* file3 = (char*) lFileName[index+2].c_str();
	  
	  lFileInterpolateCubic (file1, lFileTime[index], file2, lFileTime[index+1], file3, lFileTime[index+2], lFile, time);
	}
      else if (cntrl == 3)
	{
	  char* lFile = "Inputs/lFile";
	  char* file1 = (char*) lFileName[index-1].c_str();
	  char* file2 = (char*) lFileName[index  ].c_str();
	  char* file3 = (char*) lFileName[index+1].c_str();
	  
	  lFileInterpolateCubic (file1, lFileTime[index-1], file2, lFileTime[index], file3, lFileTime[index+1], lFile, time);
	}
    }
}

void Phase::lFileInterpolateQuadratic (char* lFile1, double time1, char* lFile2, double time2, char* lFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH], line3[MAXULFILELINELENGTH], line4[MAXULFILELINELENGTH],  line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineb[MAXULFILELINELENGTH], linec[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH], linebb[MAXULFILELINELENGTH],
       line7[MAXULFILELINELENGTH],  line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first lFile
  // ................
  FILE* file = OpenFiler (lFile1);

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

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i]) != 12)
      {
	printf ("PHASE:lFileInterpolateQuadratic: Error reading lFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second lFile
  // .................
  file = OpenFiler (lFile2);

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
  
  for (int i = 0; i < nres_2; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_2[i], &v02_2[i], &v03_2[i], &v04_2[i], &v05_2[i], &v06_2[i], &v07_2[i], &v08_2[i], &v09_2[i], &v10_2[i], &v11_2[i], &v12_2[i]) != 12)
      {
	printf ("PHASE:lFileInterpolateQuadratic: Error reading lFile_2\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate lFile data
  // ......................
  int nres;
  if (nres_2 < nres_1)
    nres = nres_2;
  else
    nres = nres_1;

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
    }

  // ........................
  // Write interpolated lFile
  // ........................
  file = OpenFilew (lFile);

  fprintf (file, "%s", line1);
  fprintf (file, "%s", line2);
  fprintf (file, "%s", line3);
  fprintf (file, "%s", line4);
  fprintf (file, "%s", line5);
  if (nres == nres_1)
    fprintf (file, "%s", lineaa);
  else
    fprintf (file, "%s", linebb);
  fprintf (file, "%s", line7);
  fprintf (file, "%s", line8);

  for (int i = 0; i < nres; i++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i]);
  
  fclose (file);

  printf ("lFile Interpolation:\n");
  printf ("%s %11.4e\n", lFile1, weight1);
  printf ("%s %11.4e\n", lFile2, weight2);
 
  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1;
  delete[] v01_2; delete[] v02_2; delete[] v03_2; delete[] v04_2;
  delete[] v05_2; delete[] v06_2; delete[] v07_2; delete[] v08_2;
  delete[] v09_2; delete[] v10_2; delete[] v11_2; delete[] v12_2;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12; 
}

void Phase::lFileInterpolateCubic (char* lFile1, double time1, char* lFile2, double time2, char* lFile3, double time3, char* lFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH], line3[MAXULFILELINELENGTH], line4[MAXULFILELINELENGTH],  line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineb[MAXULFILELINELENGTH], linec[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH], linebb[MAXULFILELINELENGTH],
       linecc[MAXULFILELINELENGTH], line7[MAXULFILELINELENGTH], line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first lFile
  // ................
  FILE* file = OpenFiler (lFile1);

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

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i]) != 12)
      {
	printf ("PHASE:lFileInterpolateCubic: Error reading lFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second lFile
  // .................
  file = OpenFiler (lFile2);

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
  
  for (int i = 0; i < nres_2; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_2[i], &v02_2[i], &v03_2[i], &v04_2[i], &v05_2[i], &v06_2[i], &v07_2[i], &v08_2[i], &v09_2[i], &v10_2[i], &v11_2[i], &v12_2[i]) != 12)
      {
	printf ("PHASE:lFileInterpolateCubic: Error reading lFile_2\n");
	exit (1);
      }

  fclose (file);

  // ................
  // Read third lFile
  // ................
  file = OpenFiler (lFile3);

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
  
  for (int i = 0; i < nres_3; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_3[i], &v02_3[i], &v03_3[i], &v04_3[i], &v05_3[i], &v06_3[i], &v07_3[i], &v08_3[i], &v09_3[i], &v10_3[i], &v11_3[i], &v12_3[i]) != 12)
      {
	printf ("PHASE:lFileInterpolateCubic: Error reading lFile_3\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate lFile data
  // ......................
  if (nres_1 == nres_2 && nres_2 != nres_3)
    {
      lFileInterpolateQuadratic (lFile1, time1, lFile2, time2, lFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2)
    {
      lFileInterpolateQuadratic (lFile2, time2, lFile3, time3, lFile, time);
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
    }

  // ........................
  // Write interpolated lFile
  // ........................
  file = OpenFilew (lFile);

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
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i]);
  
  fclose (file);

  printf ("lFile Interpolation:\n");
  printf ("%s %11.4e\n", lFile1, weight1);
  printf ("%s %11.4e\n", lFile2, weight2);
  printf ("%s %11.4e\n", lFile3, weight3);
  
  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1;
  delete[] v01_2; delete[] v02_2; delete[] v03_2; delete[] v04_2;
  delete[] v05_2; delete[] v06_2; delete[] v07_2; delete[] v08_2;
  delete[] v09_2; delete[] v10_2; delete[] v11_2; delete[] v12_2;
  delete[] v01_3; delete[] v02_3; delete[] v03_3; delete[] v04_3;
  delete[] v05_3; delete[] v06_3; delete[] v07_3; delete[] v08_3;
  delete[] v09_3; delete[] v10_3; delete[] v11_3; delete[] v12_3;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12; 
}

void Phase::lFileInterpolateQuartic (char* lFile1, double time1, char* lFile2, double time2, char* lFile3, double time3,
				     char* lFile4, double time4, char* lFile, double time)
{
  char line1[MAXULFILELINELENGTH],  line2[MAXULFILELINELENGTH],  line3[MAXULFILELINELENGTH],  line4[MAXULFILELINELENGTH], line5[MAXULFILELINELENGTH],
       linea[MAXULFILELINELENGTH],  lineb[MAXULFILELINELENGTH],  linec[MAXULFILELINELENGTH],  lined[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH],
       linebb[MAXULFILELINELENGTH], linecc[MAXULFILELINELENGTH], linedd[MAXULFILELINELENGTH], line7[MAXULFILELINELENGTH], line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first lFile
  // ................
  FILE* file = OpenFiler (lFile1);

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

  for (int i = 0; i < nres_1; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_1[i], &v02_1[i], &v03_1[i], &v04_1[i], &v05_1[i], &v06_1[i], &v07_1[i], &v08_1[i], &v09_1[i], &v10_1[i], &v11_1[i], &v12_1[i]) != 12)
      {
	printf ("PHASE::lFileInterpolateQuartic: Error reading lFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second lFile
  // .................
  file = OpenFiler (lFile2);

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
  
  for (int i = 0; i < nres_2; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_2[i], &v02_2[i], &v03_2[i], &v04_2[i], &v05_2[i], &v06_2[i], &v07_2[i], &v08_2[i], &v09_2[i], &v10_2[i], &v11_2[i], &v12_2[i]) != 12)
      {
	printf ("PHASE::lFileInterpolateQuartic: Error reading lFile_2\n");
	exit (1);
      }

  fclose (file);

  // ................
  // Read third lFile
  // ................
  file = OpenFiler (lFile3);

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
  
  for (int i = 0; i < nres_3; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_3[i], &v02_3[i], &v03_3[i], &v04_3[i], &v05_3[i], &v06_3[i], &v07_3[i], &v08_3[i], &v09_3[i], &v10_3[i], &v11_3[i], &v12_3[i]) != 12)
      {
	printf ("PHASE::lFileInterpolateQuartic: Error reading lFile_3\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read fourth lFile
  // .................
  file = OpenFiler (lFile4);

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
  
  for (int i = 0; i < nres_4; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&v01_4[i], &v02_4[i], &v03_4[i], &v04_4[i], &v05_4[i], &v06_4[i], &v07_4[i], &v08_4[i], &v09_4[i], &v10_4[i], &v11_4[i], &v12_4[i]) != 12)
      {
	printf ("PHASE::lFileInterpolateQuartic: Error reading lFile_4\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate lFile data
  // ......................
  if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 != nres_2)
    {
      lFileInterpolateQuadratic (lFile2, time2, lFile3, time3, lFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 == nres_2 && nres_4 != nres_2)
    {
      lFileInterpolateCubic (lFile1, time1, lFile2, time2, lFile3, time3, lFile, time);
      return;
    }
  else if (nres_2 == nres_3 && nres_1 != nres_2 && nres_4 == nres_2)
    {
      lFileInterpolateCubic (lFile2, time2, lFile3, time3, lFile4, time4, lFile, time);
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
    }

  // ........................
  // Write interpolated lFile
  // ........................
  file = OpenFilew (lFile);

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
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i]);
  
  fclose (file);

  printf ("lFile Interpolation:\n");
  printf ("%s %11.4e\n", lFile1, weight1);
  printf ("%s %11.4e\n", lFile2, weight2);
  printf ("%s %11.4e\n", lFile3, weight3);
  printf ("%s %11.4e\n", lFile4, weight4);
  
  // ........
  // Clean up
  // ........
  delete[] v01_1; delete[] v02_1; delete[] v03_1; delete[] v04_1;
  delete[] v05_1; delete[] v06_1; delete[] v07_1; delete[] v08_1;
  delete[] v09_1; delete[] v10_1; delete[] v11_1; delete[] v12_1;
  delete[] v01_2; delete[] v02_2; delete[] v03_2; delete[] v04_2;
  delete[] v05_2; delete[] v06_2; delete[] v07_2; delete[] v08_2;
  delete[] v09_2; delete[] v10_2; delete[] v11_2; delete[] v12_2;
  delete[] v01_3; delete[] v02_3; delete[] v03_3; delete[] v04_3;
  delete[] v05_3; delete[] v06_3; delete[] v07_3; delete[] v08_3;
  delete[] v09_3; delete[] v10_3; delete[] v11_3; delete[] v12_3;
  delete[] v01_4; delete[] v02_4; delete[] v03_4; delete[] v04_4;
  delete[] v05_4; delete[] v06_4; delete[] v07_4; delete[] v08_4;
  delete[] v09_4; delete[] v10_4; delete[] v11_4; delete[] v12_4;
  delete[] v01;   delete[] v02;   delete[] v03;   delete[] v04;
  delete[] v05;   delete[] v06;   delete[] v07;   delete[] v08;
  delete[] v09;   delete[] v10;   delete[] v11;   delete[] v12; 
}

