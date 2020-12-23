// Smoothing.cpp

#include "Flux.h"

// ###################################
// Savitsky-Gorlay smoothing algorithm
// ###################################
void Flux::Smoothing (int N, double *y)
{
  double* w = new double[N];

  for (int i = 0; i < N; i++)
    w[i] = y[i];

  double* weights = new double[5];
  for (int i = 0; i < N; i++)
    {
      int j;
      
      if (i == 0)
	{
	  weights[0] = +31.;
	  weights[1] = +9.;
	  weights[2] = -3.;
	  weights[3] = -5.;
	  weights[4] = +3.;

	  j = 2;
	}
      else if (i == 1)
	{
	  weights[0] = +9.;
	  weights[1] = +13.;
	  weights[2] = +12.;
	  weights[3] = +6.;
	  weights[4] = -5.;

	  j = 2;
	}
       else if (i == N-2)
	{
	  weights[0] = -5.;
	  weights[1] = +6.;
	  weights[2] = +12.;
	  weights[3] = +13.;
	  weights[4] = +9.;

	  j = N-3;
	}
       else if (i == N-1)
	{
	  weights[0] = +3.;
	  weights[1] = -5.;
	  weights[2] = -3.;
	  weights[3] = +9.;
	  weights[4] = +13.;

	  j = N-3;
	}
       else
	 {
	   weights[0] = -3.;
	   weights[1] = +12.;
	   weights[2] = +17.;
	   weights[3] = +12.;
	   weights[4] = -3.;

	   j = i;
	 }

      y[i] = (weights[0] * w[j-2] + weights[1] * w[j-1] + weights[2] * w[j] + weights[3] * w[j+1] + weights[4] * w[j+2]) /35.;
    }

  delete[] w;
}
