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
  double chi7[MAXSIZE];

  printf ("Reading chiph.txt\n");
  FILE* file = fopen ("chiph.txt", "r");

  int I = 1;
  while (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf",
		 &PsiN[I], &chi1[I], &chi2[I], &chi3[I],
		 &chi4[I], &chi5[I], &chi6[I], &chi7[I]) == 8)
    {
      I++;
    }

  fclose (file);

  for (int i = 1; i < I; i++)
    printf ("PsiN = %11.4e  chi_perp = (%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e)\n",
	    PsiN[i], chi1[i], chi2[i], chi3[i],
	    chi4[i], chi5[i], chi6[i], chi7[i]);

  double chi[MAXSIZE];
  for (int i = 1; i < I; i++)
    chi[i] = (chi1[i] + chi2[i] + chi3[i] + chi4[i] + chi5[i] + chi6[i] + chi7[i])/7.e4;

  PsiN[0] = 0.;
  double s0 = (PsiN[0] - PsiN[2]) * (PsiN[0] - PsiN[1]) /(PsiN[3] - PsiN[2]) /(PsiN[3] - PsiN[1]);
  double s1 = (PsiN[0] - PsiN[3]) * (PsiN[0] - PsiN[1]) /(PsiN[2] - PsiN[3]) /(PsiN[2] - PsiN[1]);
  double s2 = (PsiN[0] - PsiN[3]) * (PsiN[0] - PsiN[2]) /(PsiN[1] - PsiN[3]) /(PsiN[1] - PsiN[2]);

  chi[0] = s0 * chi[3] + s1 * chi[2] + s2 * chi[1];
  
  PsiN[I] = 1.;
  s0 = (PsiN[I] - PsiN[I-2]) * (PsiN[I] - PsiN[I-1]) /(PsiN[I-3] - PsiN[I-2]) /(PsiN[I-3] - PsiN[I-1]);
  s1 = (PsiN[I] - PsiN[I-3]) * (PsiN[I] - PsiN[I-1]) /(PsiN[I-2] - PsiN[I-3]) /(PsiN[I-2] - PsiN[I-1]);
  s2 = (PsiN[I] - PsiN[I-3]) * (PsiN[I] - PsiN[I-2]) /(PsiN[I-1] - PsiN[I-3]) /(PsiN[I-1] - PsiN[I-2]);

  chi[I] = s0 * chi[I-3] + s1 * chi[I-2] + s2 * chi[I-1];

  double chie1[MAXSIZE];
  double chie2[MAXSIZE];
  double chie3[MAXSIZE];
  double chie4[MAXSIZE];
  double chie5[MAXSIZE];
  double chie6[MAXSIZE];
  double chie7[MAXSIZE];

  printf ("Reading chi_e.txt\n");
  file = fopen ("chi_e.txt", "r");

  I = 0;
  while (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
		 &chie1[I], &chie2[I], &chie3[I],
		 &chie4[I], &chie5[I], &chie6[I], &chie7[I]) == 7)
    {
      I++;
    }

  fclose (file);

  for (int i = 0; i < I; i++)
    printf ("PsiN = %11.4e  chi_e = (%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e)\n",
	    PsiN[i], chie1[i], chie2[i], chie3[i],
	    chie4[i], chie5[i], chie6[i], chie7[i]);

  double chie[MAXSIZE];
  for (int i = 0; i < I; i++)
    {
      chie[i] = (chie1[i] + chie2[i] + chie3[i] + chie4[i] + chie5[i] + chie6[i] + chie7[i])/7.;
      if (chie[i] < 0.) chie[i] = 0.;
    }

  s0 = (PsiN[I] - PsiN[I-2]) * (PsiN[I] - PsiN[I-1]) /(PsiN[I-3] - PsiN[I-2]) /(PsiN[I-3] - PsiN[I-1]);
  s1 = (PsiN[I] - PsiN[I-3]) * (PsiN[I] - PsiN[I-1]) /(PsiN[I-2] - PsiN[I-3]) /(PsiN[I-2] - PsiN[I-1]);
  s2 = (PsiN[I] - PsiN[I-3]) * (PsiN[I] - PsiN[I-2]) /(PsiN[I-1] - PsiN[I-3]) /(PsiN[I-1] - PsiN[I-2]);

  chie[I] = s0 * chie[I-3] + s1 * chie[I-2] + s2 * chie[I-1];

  double chin1[MAXSIZE];
  double chin2[MAXSIZE];
  double chin3[MAXSIZE];
  double chin4[MAXSIZE];
  double chin5[MAXSIZE];
  double chin6[MAXSIZE];
  double chin7[MAXSIZE];

  printf ("Reading d_perp.txt\n");
  file = fopen ("d_perp.txt", "r");

  I = 0;
  while (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
		 &chin1[I], &chin2[I], &chin3[I],
		 &chin4[I], &chin5[I], &chin6[I], &chin7[I]) == 7)
    {
      I++;
    }

  fclose (file);

  for (int i = 0; i < I; i++)
    printf ("PsiN = %11.4e  D_perp = (%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e)\n",
	    PsiN[i], chin1[i], chin2[i], chin3[i],
	    chin4[i], chin5[i], chin6[i], chin7[i]);

  double chin[MAXSIZE];
  for (int i = 0; i < I; i++)
    {
      chin[i] = (chin1[i] + chin2[i] + chin3[i] + chin4[i] + chin5[i] + chin6[i] + chin7[i])/7.;
      if (chin[i] < 0.) chin[i] = 0.;
    }

  s0 = (PsiN[I] - PsiN[I-2]) * (PsiN[I] - PsiN[I-1]) /(PsiN[I-3] - PsiN[I-2]) /(PsiN[I-3] - PsiN[I-1]);
  s1 = (PsiN[I] - PsiN[I-3]) * (PsiN[I] - PsiN[I-1]) /(PsiN[I-2] - PsiN[I-3]) /(PsiN[I-2] - PsiN[I-1]);
  s2 = (PsiN[I] - PsiN[I-3]) * (PsiN[I] - PsiN[I-2]) /(PsiN[I-1] - PsiN[I-3]) /(PsiN[I-1] - PsiN[I-2]);

  chin[I] = s0 * chin[I-3] + s1 * chin[I-2] + s2 * chin[I-1];
  
  file = fopen ("c145380.2500", "w");

  fprintf (file, "%d\n", I+1);
  for (int i = 0; i <= I; i++)
    fprintf (file, "%19.6e  %19.6e  %19.6e  %19.6e\n", PsiN[i], chi[i], chie[i], chin[i]);

  fclose (file);
 }
