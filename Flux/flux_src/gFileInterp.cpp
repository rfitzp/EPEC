// gFileInterpolate.cpp

#include "Flux.h"

// ###############################
// Functions to interpolate gFiles
// ###############################
void Flux::_gFileInterpCubic (char* gFile1, double time1, char* gFile2, double time2, char* gFile3, double time3,
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

 gFileInterpolateCubic ();
}

void Flux::_gFileInterpQuartic (char* gFile1, double time1, char* gFile2, double time2, char* gFile3, double time3,
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

 gFileInterpolateQuartic ();
}

void Flux::gFileInterp (vector<string> gFileName, vector<double> gFileTime, int gFileNumber, double time)
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
      char* gFile = "gFile";
      char* file1 = (char*) gFileName[index-1].c_str();
      char* file2 = (char*) gFileName[index  ].c_str();
      char* file3 = (char*) gFileName[index+1].c_str();
      char* file4 = (char*) gFileName[index+2].c_str();
      
      _gFileInterpQuartic (file1, gFileTime[index-1], file2, gFileTime[index], file3, gFileTime[index+1],
			   file4, gFileTime[index+2], gFile, time);
    }
  else if (cntrl == 2)
    {
      char* gFile = "gFile";
      char* file1 = (char*) gFileName[index  ].c_str();
      char* file2 = (char*) gFileName[index+1].c_str();
      char* file3 = (char*) gFileName[index+2].c_str();
      
      _gFileInterpCubic (file1, gFileTime[index], file2, gFileTime[index+1], file3, gFileTime[index+2], gFile, time);
    }
  else if (cntrl == 3)
    {
      char* gFile = "gFile";
      char* file1 = (char*) gFileName[index-1].c_str();
      char* file2 = (char*) gFileName[index  ].c_str();
      char* file3 = (char*) gFileName[index+1].c_str();
      
      _gFileInterpCubic (file1, gFileTime[index-1], file2, gFileTime[index], file3, gFileTime[index+1], gFile, time);
    }
}
