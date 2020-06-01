import graph;
     
size(500,500,IgnoreAspect);

file    in = input("rkminusne.txt").line();
real[][] A = in.dimension (0,0);
int nres   = A[0].length - 1;
A          = transpose(A);
int ntim   = A[0].length;

file    in1 = input("rkplusne.txt").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

int nav    = 100;
     
real[] p  = A[0];
real[] p1 = A1[0];
pen s;

for (int j = 0; j < ntim; ++j)
  {
    if (j > nav)
      {
	for (int i = 0; i < nres; ++i)
	  {
	    real sum = 0.; real sum1 = 0.;
	    for (int k = 0; k < nav; ++k)
	      {
		sum  += A[i+1][j-k];
		sum1 += A1[i+1][j-k];	
	      }
	    sum  /= (int) nav;
	    sum1 /= (int) nav;
	    A[i+1][j]  = sum;
	    A1[i+1][j] = sum1;
	  }
      }
    else
      {
	for (int i = 0; i < nres; ++i)
	  {
	    real sum = 0.; real sum1 = 0.;
	    for (int k = 0; k < nav; ++k)
	      {
		sum  += A[i+1][k];
		sum1 += A1[i+1][k];
	      }
	    sum /= (int) nav;
	    sum1 /= (int) nav;
	    A[i+1][j]  = sum;
	    A1[i+1][j] = sum1;
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
 
    draw (graph (p,  A[i+1]),  s);
    draw (graph (p1, A1[i+1]), s);
  }
s = dotted + black + 1;
ylimits (0.4, 1.0, Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$r/a$",LeftRight,RightTicks);
