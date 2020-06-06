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

      char           Basename[MAXFILENAMELENGTH];
      char           Filename[MAXFILENAMELENGTH];
      char           filename[MAXFILENAMELENGTH];
      vector<string> gFileName;
      double         filetime;
      vector<double> gFileTime;
      int            gFileNumber;
  
      printf ("Reading gFile data:\n");

      FILE* file = OpenFiler ((char*) "gFileIndex");

      fscanf (file, "%s", &Basename);

      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  strcpy (Filename, Basename);
	  strcat (Filename, filename);
	  
	  gFileName.push_back (Filename);
	  gFileTime.push_back (filetime);
	}
      gFileNumber = gFileTime.size ();

      fclose (file);
      
      gFileInterp (gFileName, gFileTime, gFileNumber, TIME);
    }

  // Read gFile
  gFileRead ();
}

