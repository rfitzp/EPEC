// Field.cpp

#include "Field.h"

Field::Field ()
{
  N = 0;
}

Field::Field (int n)
{
  N    = n;
  X    = new double[N];
  Y    = new double[N];
  dYdX = new double[N];
}

Field::~Field ()
{
  if (N > 0)
    {
      delete[] X;
      delete[] Y;
      delete[] dYdX;
    }
}

void Field::resize (int n)
{
  N    = n;
  X    = new double[N];
  Y    = new double[N];
  dYdX = new double[N];
}

int Field::GetN ()
{
  return N;
}

void Field::PushData (int i, double x, double y, double dydx)
{
  X[i]    = x;
  Y[i]    = y;
  dYdX[i] = dydx;
}

void Field::PullData (int i, double& x, double& y, double& dydx)
{
  x    = X[i];
  y    = Y[i];
  dydx = dYdX[i];
}

double Field::GetX (int i)
{
  return X[i];
}

double Field::GetY (int i)
{
  return Y[i];
}

double Field::GetdYdX (int i)
{
  return dYdX[i];
}



