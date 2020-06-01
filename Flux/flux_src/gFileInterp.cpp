// gFileInterpolate.cpp

#include "Flux.h"

// ###############################
// Functions to interpolate gFiles
// ###############################
void Flux::_gFileInterp (char* gFile1, double time1, char* gFile2, double time2, char* gFile3, double time3, char* gFile, double time)
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

 gFileInterpolate ();
}

void Flux::gFileInterp (vector<string> gFileName, vector<double> gFileTime, int gFileNumber, double time)
{
  int    index;
  double _time;

  if (time < gFileTime[0])
    {
      index = 0;
      _time = gFileTime[0];
    }
  else if (time >= gFileTime[gFileNumber-1])
    {
      index = gFileNumber - 3;
      _time = gFileTime[gFileNumber-1];
    }
  else
    {
      for (int i = 0; i < gFileNumber-1; i++)
	if (time >= gFileTime[i] && time < gFileTime[i+1])
	  {
	    index = i;
	    _time = time;

	    if (index > gFileNumber-3)
	      index = gFileNumber - 3;
	  }
    }

  char* gFile = "gFile";
  char* file1 = (char*) gFileName[index  ].c_str();
  char* file2 = (char*) gFileName[index+1].c_str();
  char* file3 = (char*) gFileName[index+2].c_str();
  
  _gFileInterp (file1, gFileTime[index], file2, gFileTime[index+1], file3, gFileTime[index+2], gFile, _time);
}
