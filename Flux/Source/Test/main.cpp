#include <stdio.h>

double ToroidalP (int m, int n, double z);

int main ()
{
  int m, n;
  double z;

  printf ("m n z ?? ");
  scanf ("%d %d %lf", &m, &n, &z);

  printf ("P = %17.10e\n", ToroidalP (m, n, z));
}


