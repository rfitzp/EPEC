// Smoothing.cpp

#include "Flux.h"

// ###################################
// Savitsky-Gorlay smoothing algorithm
// ###################################
void Flux::Smoothing (int N, double *y)
{
  double* w       = new double[N];
  double* weights = new double[5];

  for (int i = 0; i < N; i++)
    w[i] = y[i];

  for (int i = 2; i < N-2; i++)
    {
      weights[0] = -3.;
      weights[1] = +12.;
      weights[2] = +17.;
      weights[3] = +12.;
      weights[4] = -3.;

      y[i] = (weights[0] * w[i-2] + weights[1] * w[i-1] + weights[2] * w[i] + weights[3] * w[i+1] + weights[4] * w[i+2]) /35.;
    }

  delete[] w;
  delete[] weights;
}
