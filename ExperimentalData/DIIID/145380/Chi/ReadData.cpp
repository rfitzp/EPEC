#include <stdio.h>

#define MAXSIZE 100

int main ()
{
  double PsiN[MAXSIZE];
  double chi1[MAXSIZE];
  double chi2[MAXSIZE];
  double chi3[MAXSIZE];
  double chi4[MAXSIZE];
  double chi5[MAXSIZE];
  double chi6[MAXSIZE];

  printf ("Reading coeff.txt\n");
  FILE* file = fopen ("coeff.txt", "r");

  int I = 1;
  while (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
		 &PsiN[I], &chi1[I], &chi2[I], &chi3[I],
		 &chi4[I], &chi5[I], &chi6[I]) == 7)
    {
      I++;
    }

  fclose (file);

  PsiN[0] = 0.;
  double s0 = (PsiN[0] - PsiN[2]) * (PsiN[0] - PsiN[1]) /(PsiN[3] - PsiN[2]) /(PsiN[3] - PsiN[1]);
  double s1 = (PsiN[0] - PsiN[3]) * (PsiN[0] - PsiN[1]) /(PsiN[2] - PsiN[3]) /(PsiN[2] - PsiN[1]);
  double s2 = (PsiN[0] - PsiN[3]) * (PsiN[0] - PsiN[2]) /(PsiN[1] - PsiN[3]) /(PsiN[1] - PsiN[2]);

  chi1[0] = s0 * chi1[3] + s1 * chi1[2] + s2 * chi1[1];
  chi3[0] = s0 * chi3[3] + s1 * chi3[2] + s2 * chi3[1];
  chi5[0] = s0 * chi5[3] + s1 * chi5[2] + s2 * chi5[1];
     
  for (int i = 0; i < I; i++)
    printf ("PsiN = %11.4e  D_perp = %11.4e  chi_perp = %11.4e  chi_phi = %11.4e\n",
	    PsiN[i], chi1[i], chi3[i], chi5[i]);

  file = fopen ("c145380.2500", "w");

  fprintf (file, "%d\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%19.6e  %19.6e  %19.6e  %19.6e\n", PsiN[i], chi5[i], chi3[i], chi1[i]);

  fclose (file);
 }
