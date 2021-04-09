// gFileInterpolate.cpp

// PROGRAM ORGANIZATION:
//
// void Flux:: gFileInterp               (vector<string> gFileName,   vector<double> gFileTime,   int gFileNumber,  double time)
// void Flux:: gFileInterpolateLinearv   (char* gFile1, double time1, char* gFile,  double time)
// void Flux:: gFileInterpolateQuadratic (char* gFile1, double time1, char* gFile2, double time2, char* gFile,      double time)
// void Flux:: gFileInterpolateCubic     (char* gFile1, double time1, char* gFile2, double time2, char* gFile3,     double time3, char* gFile, double time)
// void Flux:: gFileInterpolateQuartic   (char* gFile1, double time1, char* gFile2, double time2, char* gFile3,     double time3,
//				          char* gFile4, double time4, char* gFile,  double time)

#include "Flux.h"

// ###############################
// Functions to interpolate gFiles
// ###############################
void Flux::gFileInterp (vector<string> gFileName, vector<double> gFileTime, int gFileNumber, double time)
{
  if (gFileNumber < 1)
    {
      printf ("FLUX::gFileInterp - gFileNumber must be greater than zero\n");
      exit (1);
    }
  else if (gFileNumber == 1)
    {
      char* gFile = "Inputs/gFile";
      char* file1 = (char*) gFileName[0].c_str();

      gFileInterpLinear (file1, gFileTime[0], gFile, time);
    }
  else if (gFileNumber == 2)
    {
      char* gFile = "Inputs/gFile";
      char* file1 = (char*) gFileName[0].c_str();
      char* file2 = (char*) gFileName[1].c_str();

      gFileInterpQuadratic (file1, gFileTime[0], file2, gFileTime[1], gFile, time);
    }
  else if (gFileNumber == 3)
    {
      char* gFile = "Inputs/gFile";
      char* file1 = (char*) gFileName[0].c_str();
      char* file2 = (char*) gFileName[1].c_str();
      char* file3 = (char*) gFileName[2].c_str();

      gFileInterpCubic (file1, gFileTime[0], file2, gFileTime[1], file3, gFileTime[2], gFile, time);
    }
  else if (gFileNumber == 4)
    {
      char* gFile = "Inputs/gFile";
      char* file1 = (char*) gFileName[0].c_str();
      char* file2 = (char*) gFileName[1].c_str();
      char* file3 = (char*) gFileName[2].c_str();
      char* file4 = (char*) gFileName[3].c_str();

      gFileInterpQuartic (file1, gFileTime[0], file2, gFileTime[1], file3, gFileTime[2],
			  file4, gFileTime[3], gFile, time);
    }
  else
    {
      int index, cntrl;
      
      if (time < gFileTime[0])
	{
	  index = 0;
	  cntrl = 2;
	}
      else if (time >= gFileTime[gFileNumber-1])
	{
	  index = gFileNumber - 2;
	  cntrl = 3;
	}
      else
	{
	  for (int i = 0; i < gFileNumber-1; i++)
	    if (time >= gFileTime[i] && time < gFileTime[i+1])
	      {
		index = i;
		
		if (index == 0)
		  cntrl = 2;
		else if (index == gFileNumber-2)
		  cntrl = 3;
		else
		  cntrl = 1;
	      }
	}
      
      if (cntrl == 1)
	{
	  char* gFile = "Inputs/gFile";
	  char* file1 = (char*) gFileName[index-1].c_str();
	  char* file2 = (char*) gFileName[index  ].c_str();
	  char* file3 = (char*) gFileName[index+1].c_str();
	  char* file4 = (char*) gFileName[index+2].c_str();
	  
	  gFileInterpQuartic (file1, gFileTime[index-1], file2, gFileTime[index], file3, gFileTime[index+1],
			      file4, gFileTime[index+2], gFile, time);
	}
      else if (cntrl == 2)
	{
	  char* gFile = "Inputs/gFile";
	  char* file1 = (char*) gFileName[index  ].c_str();
	  char* file2 = (char*) gFileName[index+1].c_str();
	  char* file3 = (char*) gFileName[index+2].c_str();
	  
	  gFileInterpCubic (file1, gFileTime[index], file2, gFileTime[index+1], file3, gFileTime[index+2], gFile, time);
	}
      else if (cntrl == 3)
	{
	  char* gFile = "Inputs/gFile";
	  char* file1 = (char*) gFileName[index-1].c_str();
	  char* file2 = (char*) gFileName[index  ].c_str();
	  char* file3 = (char*) gFileName[index+1].c_str();
	  
	  gFileInterpCubic (file1, gFileTime[index-1], file2, gFileTime[index], file3, gFileTime[index+1], gFile, time);
	}
    }
}

void Flux::gFileInterpLinear (char* gFile1, double time1, char* gFile, double time)
{
  FILE* file = OpenFilew ("Interface.txt");

  fprintf (file, "%s\n",     gFile1);
  fprintf (file, "%16.9e\n", time1);
  fprintf (file, "%s\n",     gFile);
  fprintf (file, "%16.9e\n", time);

  fclose (file);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  double weight1 = 1.;
  fprintf (monitor, "gFile Interpolation:\n");
  fprintf (monitor, "%s %16.9e\n", gFile1, weight1);
  fclose (monitor);
  
  gFileInterpolateLinear ();
}


void Flux::gFileInterpQuadratic (char* gFile1, double time1, char* gFile2, double time2, char* gFile, double time)
{
  FILE* file = OpenFilew ("Interface.txt");

  fprintf (file, "%s\n",     gFile1);
  fprintf (file, "%16.9e\n", time1);
  fprintf (file, "%s\n",     gFile2);
  fprintf (file, "%16.9e\n", time2);
  fprintf (file, "%s\n",     gFile);
  fprintf (file, "%16.9e\n", time);

  fclose (file);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  double weight1 = (time - time2) /(time1 - time2);
  double weight2 = (time - time1) /(time2 - time1);
  fprintf (monitor, "gFile Interpolation:\n");
  fprintf (monitor, "%s %16.9e\n", gFile1, weight1);
  fprintf (monitor, "%s %16.9e\n", gFile2, weight2);
  fclose (monitor);
  
  gFileInterpolateQuadratic ();
}

void Flux::gFileInterpCubic (char* gFile1, double time1, char* gFile2, double time2, char* gFile3, double time3,
			     char* gFile, double time)
{
  FILE* file = OpenFilew ("Interface.txt");

  fprintf (file, "%s\n",     gFile1);
  fprintf (file, "%16.9e\n", time1);
  fprintf (file, "%s\n",     gFile2);
  fprintf (file, "%16.9e\n", time2);
  fprintf (file, "%s\n",     gFile3);
  fprintf (file, "%16.9e\n", time3);
  fprintf (file, "%s\n",     gFile);
  fprintf (file, "%16.9e\n", time);

  fclose (file);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  double weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3);
  double weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3);
  double weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2);
  fprintf (monitor, "gFile Interpolation:\n");
  fprintf (monitor, "%s %16.9e\n", gFile1, weight1);
  fprintf (monitor, "%s %16.9e\n", gFile2, weight2);
  fprintf (monitor, "%s %16.9e\n", gFile3, weight3);
  fclose (monitor);

  gFileInterpolateCubic ();
}

void Flux::gFileInterpQuartic (char* gFile1, double time1, char* gFile2, double time2, char* gFile3, double time3,
			       char* gFile4, double time4, char* gFile, double time)
{
  FILE* file = OpenFilew ("Interface.txt");

  fprintf (file, "%s\n",     gFile1);
  fprintf (file, "%16.9e\n", time1);
  fprintf (file, "%s\n",     gFile2);
  fprintf (file, "%16.9e\n", time2);
  fprintf (file, "%s\n",     gFile3);
  fprintf (file, "%16.9e\n", time3);
  fprintf (file, "%s\n",     gFile4);
  fprintf (file, "%16.9e\n", time4);
  fprintf (file, "%s\n",     gFile);
  fprintf (file, "%16.9e\n", time);

  fclose (file);

  FILE* monitor = OpenFilea ((char*) "../IslandDynamics/Outputs/monitor.txt");
  double weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4);
  double weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4);
  double weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4);
  double weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3);
  fprintf (monitor, "gFile Interpolation:\n");
  fprintf (monitor, "%s %16.9e\n", gFile1, weight1);
  fprintf (monitor, "%s %16.9e\n", gFile2, weight2);
  fprintf (monitor, "%s %16.9e\n", gFile3, weight3);
  fprintf (monitor, "%s %16.9e\n", gFile4, weight4);
  fclose (monitor);

  gFileInterpolateQuartic ();
}

