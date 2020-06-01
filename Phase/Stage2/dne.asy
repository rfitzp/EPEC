import graph;
     
size(500,500,IgnoreAspect);

file    in = input("dne.txt").line();
real[][] A = in.dimension (0,0);
int nres   = A[0].length - 1;
A          = transpose(A);
int ntim   = A[0].length;

int nav    = 100;
     
real[] p = A[0];
pen s;

for (int j = 0; j < ntim; ++j)
  {
    if (j > nav)
      {
	for (int i = 0; i < nres; ++i)
	  {
	    real sum = 0.;
	    for (int k = 0; k < nav; ++k)
	      {
		sum += A[i+1][j-k];
	      }
	    sum /= (int) nav;
	    A[i+1][j] = sum;
	  }
      }
    else
      {
	for (int i = 0; i < nres; ++i)
	  {
	    real sum = 0.;
	    for (int k = 0; k < nav; ++k)
	      {
		sum += A[i+1][k];
	      }
	    sum /= (int) nav;
	    A[i+1][j] = sum;
	  }
      }
  }

for (int i = 0; i < nres; i += 1)
  {
    if (i%4 == 0)
      s = black;
      else if (i%4 == 1)
	s = red;
      else if (i%4 == 2)
	s = green;
      else if (i%4 == 3)
	s = blue;
 
    draw (graph (p, A[i+1]), s);
  }
s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\delta n_e (10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
