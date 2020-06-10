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
      // Save pwd
      char pwd[MAXFILENAMELENGTH];
      getcwd (pwd, MAXFILENAMELENGTH);

      // Remove gFile
      chdir ("Inputs");
      system ("rm -rf gFile");

      // Get gFiles directory
      char gFileDir[MAXFILENAMELENGTH];
      readlink ("gFiles", gFileDir, MAXFILENAMELENGTH);
      chdir (pwd);

      // Read gFile data
      char           Basename[MAXFILENAMELENGTH];
      char           Filename[MAXFILENAMELENGTH];
      char           filename[MAXFILENAMELENGTH];
      vector<string> gFileName;
      double         filetime;
      vector<double> gFileTime;
      int            gFileNumber;
  
      printf ("Reading gFile data:\n");

      chdir (gFileDir);
      getcwd (Basename, MAXFILENAMELENGTH);
      strcat (Basename, "/");

      FILE* file = OpenFiler ((char*) "Index");

      while (fscanf (file, "%s %lf", &filename, &filetime) == 2)
	{
	  strcpy (Filename, Basename);
	  strcat (Filename, filename);
	  
	  gFileName.push_back (Filename);
	  gFileTime.push_back (filetime);
	}
      gFileNumber = gFileTime.size ();

      fclose (file);
      chdir (pwd);

      // Interpolate gFiles
      gFileInterp (gFileName, gFileTime, gFileNumber, TIME);
    }
 
  // Read gFile
  gFileRead ();
}

