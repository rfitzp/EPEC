#include <stdio.h>
#include <stdlib.h>

int main ()
{
  double Psin[500], q[500], Te[500], ne[500], Ti[500], ni[500];
  double c; int d;
  char   line1[1000];

  FILE*  file = fopen ("18594_6450ms_profiles.csv", "r");
 
  fgets (line1, 1000, file);

  for (int i = 0; i < 402; i++)
    {
      fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
      	      &Psin[i], &c, &c, &c, &c, &c, &q[i], &Te[i], &ne[i], &Ti[i], &ni[i], &c, &c, &c, &c, &d, &d, &d);
      printf ("PsiN = %11.4e q = %11.4e Te = %11.4e ne = %11.4e  Ti = %11.4e  ni = %11.4e\n", Psin[i], q[i], Te[i], ne[i], Ti[i], ni[i]);
    }

  fclose (file);

  double Zeff = 1.833333333, Z = 6.;

  double nI[500], nb[500];

  for (int i = 0; i < 402; i++)
    {
      nI[i] = (Zeff * ne[i] - ni[i]) /Z/Z;
      nb[i] = ne[i] - ni[i] - Z * nI[i];
      
      printf ("PsiN = %11.4e nI/ni = %11.4e nb/ni = %11.4e\n", Psin[i], nI[i]/ni[i], nb[i]/ni[i]);
    }
}
