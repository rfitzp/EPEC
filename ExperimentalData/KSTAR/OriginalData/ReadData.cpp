#include <stdio.h>
#include <stdlib.h>

int main ()
{
  double a1;

  printf ("a1 ?? ");
  scanf  ("%lf", &a1);
  printf ("a1 = %11.4e\n", a1);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double Psin[500], R[500], Te[500], ne[500], Ti[500], Vtor[500];
  char   line1[1000];
  int    d, I = 401;
  double c;
  FILE*  file = fopen ("18594_6450ms_profiles.txt", "r");
 
  fgets (line1, 1000, file);

  for (int i = 0; i < I; i++)
    {
      fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
      	      &Psin[i], &c, &R[i], &c, &c, &c, &c, &Te[i], &ne[i], &Ti[i], &c, &Vtor[i], &c, &c, &c, &d, &d, &d);

      ne[i] /= 10.;
      ne[i] *= a1;
      Te[i] *= a1;
      Ti[i] *= a1;
 
      printf ("PsiN = %11.4e R = %11.4e Te = %11.4e ne = %11.4e Ti = %11.4e Vtor = %11.4e\n", Psin[i], R[i], Te[i], ne[i], Ti[i], Vtor[i]);
    }

  fclose (file);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double ni[500], nI[500], nb[500], wtor[500];
  double Zeff = 2., ZI = 6.;

  for (int i = 0; i < I; i++)
    {
      ni[i] = ((ZI - Zeff) /(ZI - 1.))     * ne[i];
      nI[i] = ((Zeff - 1.) /(ZI - 1.) /ZI) * ne[i];
      nb[i] = 0.;
      
      wtor[i] = Vtor[i] /R[i];
    }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double dTe[500], dne[500], dTi[500], dni[500], dnI[500], dwtor[500];
  
  for (int i = 1; i < I-1; i++)
    {
      dne[i]   = (ne  [i+1] - ne  [i-1]) /(Psin[i+1] - Psin[i-1]);
      dTe[i]   = (Te  [i+1] - Te  [i-1]) /(Psin[i+1] - Psin[i-1]);
      dni[i]   = (ni  [i+1] - ni  [i-1]) /(Psin[i+1] - Psin[i-1]);
      dTi[i]   = (Ti  [i+1] - Ti  [i-1]) /(Psin[i+1] - Psin[i-1]);
      dnI[i]   = (nI  [i+1] - nI  [i-1]) /(Psin[i+1] - Psin[i-1]);
      dwtor[i] = (wtor[i+1] - wtor[i-1]) /(Psin[i+1] - Psin[i-1]);
    }
  dne  [0]   = 2.*dne  [1]   - dne  [2];
  dTe  [0]   = 2.*dTe  [1]   - dTe  [2];
  dni  [0]   = 2.*dni  [1]   - dni  [2];
  dTi  [0]   = 2.*dTi  [1]   - dTi  [2];
  dnI  [0]   = 2.*dnI  [1]   - dnI  [2];
  dwtor[0]   = 2.*dwtor[1]   - dwtor[2];
  dne  [I-1] = 2.*dne  [I-2] - dne  [I-3];
  dTe  [I-1] = 2.*dTe  [I-2] - dTe  [I-3];
  dni  [I-1] = 2.*dni  [I-2] - dni  [I-3];
  dTi  [I-1] = 2.*dTi  [I-2] - dTi  [I-3];
  dnI  [I-1] = 2.*dnI  [I-2] - dnI  [I-3];
  dwtor[I-1] = 2.*dwtor[I-2] - dwtor[I-3];

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double psin[2000], web[2000], dweb[2000];
  int    J = 1533;

  file = fopen ("18594_6450ms_web.txt", "r");
 
  fgets (line1, 1000, file);

  for (int i = 0; i < J; i++)
    {
      fscanf (file, "%lf %lf", &psin[i], &web[i]);
      web[i] /= 1.e3;
    }

  fclose (file);

  for (int i = 1; i < J-1; i++)
    {
      dweb[i] = (web[i+1] - web[i-1]) /(psin[i+1] - psin[i-1]);
    }
  dweb[0]   = 2.*dweb[1]   - dweb[2];
  dweb[J-1] = 2.*dweb[J-2] - dweb[J-3];

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file = fopen ("p18594.6450", "w");

  fprintf (file, "%3d psinorm ne(10^20/m^3) dnedpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], ne[i], dne[i]);

  fprintf (file, "%3d psinorm te(KeV) dtedpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], Te[i], dTe[i]);

  fprintf (file, "%3d psinorm ni(10^20/m^3) dnidpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], ni[i], dni[i]);

  fprintf (file, "%3d psinorm ti(KeV) dtidpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], Ti[i], dTi[i]);

  fprintf (file, "%3d psinorm nb(10^20/m^3) dnbdpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], nb[i], nb[i]);

  fprintf (file, "%3d psinorm omeg(kRad/s) domeg/dpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], wtor[i], dwtor[i]);

  fprintf (file, "%3d psinorm omgeb(kRad/s) domgeb/dpsiN\n", J);
  for (int i = 0; i < J; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", psin[i], web[i], dweb[i]);

  fprintf (file, "%3d psinorm nz1(10^20/m^3) dnz1dpsiN\n", I);
  for (int i = 0; i < I; i++)
    fprintf (file, "%10.6f %10.6f %10.6f\n", Psin[i], nI[i], dnI[i]);

  int    N = 3;
  double Z = 1., A = 2., AI = 12.010700;

  fprintf (file, "%3d N Z A of ION SPECIES\n", N);
  fprintf (file, "%10.6f %10.6f %10.6f\n",     ZI, ZI, AI);
  fprintf (file, "%10.6f %10.6f %10.6f\n",     Z,  Z,  A);
  fprintf (file, "%10.6f %10.6f %10.6f\n",     Z,  Z,  A);
  
  fclose (file);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double ln[500], le[500], lp[500];
  
  for (int i = 1; i < I-1; i++)
    {
      lp[i] = 1. /(Vtor[i+1] - Vtor[i]);
      le[i] = 1. /(Te[i+1]   - Te[i]);
      ln[i] = 1. /(ne[i+1]   - ne[i]);
    }
  lp[0]   = 2.*lp[1]   - lp[2];
  le[0]   = 2.*le[1]   - le[2];
  ln[0]   = 2.*ln[1]   - ln[2];
  lp[I-1] = 2.*lp[I-2] - lp[I-3];
  le[I-1] = 2.*le[I-2] - le[I-3];
  ln[I-1] = 2.*ln[I-2] - ln[I-3];
 
  double chi = 1.;
  
  file = fopen ("c18594.6450", "w");

  fprintf (file, "%d\n", I);
  for (int i = 0; i < I; i++)
    {
      fprintf (file, "%18.6e %18.6e %18.6e %18.6e\n", Psin[i], chi*lp[0]/lp[i], chi*le[0]/le[i], chi*ln[0]/ln[i]);
    }

  fclose (file);
}
