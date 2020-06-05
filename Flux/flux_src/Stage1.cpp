// Stage1.cpp

#include "Flux.h"

// ###################################################
// Function to input gFile data and output Stage1 data
// ###################################################
void Flux::Stage1 ()
{
  // Interpolate gFiles
  if (INTG > 0 && TIME > 0.)
    {
      system ("rm -rf gFile");

      char           filename[MAXFILENAMELENGTH];
      vector<string> gFileName;
      double         filetime;
      vector<double> gFileTime;
      int            gFileNumber = 0;
      
      printf ("Reading gFile data:\n");

      FILE* file = OpenFiler ((char*) "gFileIndex");

      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  gFileName.push_back (filename);
	  gFileTime.push_back (filetime);
	  gFileNumber++;
	}
      gFileNumber--;

      fclose (file);
      
      gFileInterp (gFileName, gFileTime, gFileNumber, TIME);
    }

  // Read gFile
  gFileRead ();
}

