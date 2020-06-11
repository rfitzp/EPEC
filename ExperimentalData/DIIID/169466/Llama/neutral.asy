import graph;

real Rmin = 1.8;
real Rmax = 1.96;
real Rsep = 1.92;
real nmax = 1.4;
real len  = 0.015;

real f (real R)
{
return nmax / (1. + (R - Rsep) * (R - Rsep) /len/len);
}

size (500, 500, IgnoreAspect);

file    in = input ("datapoints").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] x  = A[0];
real[] y  = A[1];

pen s = solid + 1.5;
draw (graph (x, y), s, marker(scale(1.5mm)*polygon(4)));

s = dotted + red + 1.5;

draw (graph (f, Rmin, Rmax), s); 

s = dotted + 1.5;
limits ((1.8, 0.), (1.96, 1.5), Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("$R({\rm m})$", BottomTop, LeftTicks);
yaxis("$n_n(10^{17}\,{\rm m}^{-3})$",      LeftRight,  RightTicks);
