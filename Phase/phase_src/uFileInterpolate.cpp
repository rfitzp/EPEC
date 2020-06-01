// uFileInterpolate.h

#include "Phase.h"

// ###############################
// Functions to interpolate uFiles
// ###############################
void Phase::uFileInterpolate (char* uFile1, double time1, char* uFile2, double time2, char* uFile3, double time3, char* uFile, double time)
{
  char line1[MAXULFILELINELENGTH], line2[MAXULFILELINELENGTH], line3[MAXULFILELINELENGTH], line4[MAXULFILELINELENGTH], line5[MAXULFILELINELENGTH],
    linea[MAXULFILELINELENGTH], lineb[MAXULFILELINELENGTH], linec[MAXULFILELINELENGTH], lineaa[MAXULFILELINELENGTH], linebb[MAXULFILELINELENGTH],
    linecc[MAXULFILELINELENGTH], line7[MAXULFILELINELENGTH], line8[MAXULFILELINELENGTH];
  
  // ................
  // Read first uFile
  // ................
  FILE* file = OpenFiler (uFile1);

  fgets (line1, MAXULFILELINELENGTH, file);
  fgets (line2, MAXULFILELINELENGTH, file);
  fgets (line3, MAXULFILELINELENGTH, file);
  fgets (line4, MAXULFILELINELENGTH, file);
  fgets (line5, MAXULFILELINELENGTH, file);
  fgets (linea, MAXULFILELINELENGTH, file);
  fgets (line7, MAXULFILELINELENGTH, file);
  fgets (line8, MAXULFILELINELENGTH, file);
  strcpy (lineaa, linea);

  char* token = strtok (linea, " "); 
  token = strtok (NULL, " "); 
  token = strtok (NULL, " ");
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
	printf ("Error reading uFile_1\n");
	exit (1);
      }

  fclose (file);

  // .................
  // Read second uFile
  // .................
  file = OpenFiler (uFile2);

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
  token = strtok (NULL, " "); 
  token = strtok (NULL, " ");
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
	printf ("Error reading uFile_2\n");
	exit (1);
      }

  fclose (file);

  // ................
  // Read third uFile
  // ................
  file = OpenFiler (uFile3);

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
  token = strtok (NULL, " "); 
  token = strtok (NULL, " ");
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
	printf ("Error reading uFile_3\n");
	exit (1);
      }

  fclose (file);

  // ......................
  // Interpolate uFile data
  // ......................
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
  // Write interpolated uFile
  // ........................
  file = OpenFilew (uFile);

  fprintf (file, "%s", line1);
  fprintf (file, "%s", line2);
  fprintf (file, "%s", line3);
  fprintf (file, "%s", line4);
  fprintf (file, "%s", line5);
  if (nres_1 <= nres_2 && nres_1 <= nres_3)
    fprintf (file, "%s", lineaa);
  else if (nres_2 <= nres_1 && nres_2 <= nres_3)
    fprintf (file, "%s", linebb);
  else
    fprintf (file, "%s", linecc);
  fprintf (file, "%s", line7);
  fprintf (file, "%s", line8);

  for (int i = 0; i < nres; i++)
    fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     v01[i], v02[i], v03[i], v04[i], v05[i], v06[i], v07[i], v08[i], v09[i], v10[i], v11[i], v12[i]);
  
  fclose (file);

  printf ("uFile Interpolation:\n");
  printf ("%s %11.4e\n", uFile1, weight1);
  printf ("%s %11.4e\n", uFile2, weight2);
  printf ("%s %11.4e\n", uFile3, weight3);
  
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

void Phase::uFileInterp (vector<string> uFileName, vector<double> uFileTime, vector<int> uFileNtor, vector<int> uFileMmin, vector<int> uFileMmax, int uFileNumber,
			 double time, int &Ntor, int& Mmin, int& Mmax)
{
  int    index;
  double _time;

  if (time < uFileTime[0])
    {
      index = 0;
      _time = uFileTime[0];
      Ntor  = uFileNtor[0];
      Mmin  = uFileMmin[0];
      Mmax  = uFileMmax[0];
    }
  else if (time >= uFileTime[uFileNumber-1])
    {
      index = uFileNumber - 3;
      _time = uFileTime[uFileNumber-1];
      Ntor  = uFileNtor[uFileNumber-1];
      Mmin  = uFileMmin[uFileNumber-1];
      Mmax  = uFileMmax[uFileNumber-1];
    }
  else
    {
      for (int i = 0; i < uFileNumber-1; i++)
	if (time >= uFileTime[i] && time < uFileTime[i+1])
	  {
	    index = i;
	    _time = time;

	    if (index > uFileNumber-3)
	      index = uFileNumber - 3;

	    if (uFileMmax[i+1] > uFileMmax[i])
	      {
		Ntor = uFileNtor[i];
		Mmin = uFileMmin[i];
		Mmax = uFileMmax[i];
	      }
	    else
	      {
		Ntor = uFileNtor[i+1];
		Mmin = uFileMmin[i+1];
		Mmax = uFileMmax[i+1];
	      }
	    if (uFileMmax[i+2] < Mmax)
	      {
		Ntor = uFileNtor[i+2];
		Mmin = uFileMmin[i+2];
		Mmax = uFileMmax[i+2];
	      }
	  }
    }

  char* uFile = "uFile";
  char* file1 = (char*) uFileName[index  ].c_str();
  char* file2 = (char*) uFileName[index+1].c_str();
  char* file3 = (char*) uFileName[index+2].c_str();
  
  uFileInterpolate (file1, uFileTime[index], file2, uFileTime[index+1], file3, uFileTime[index+2], uFile, _time);
}

