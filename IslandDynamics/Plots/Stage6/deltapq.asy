import graph;

size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/deltap.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

real[] t = A[0];
real[] q = A[1];
real[] p = A[2];

pen s = black + dotted + 0.5;

draw (graph (q, p), s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, black));

s = black + 1 + dotted;

yequals (0., s);

ylimits (-0.2,0.,Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$q_{95}$",                BottomTop, RightTicks(n=5));
yaxis ("$-{\mit\Delta}P/P_{\rm ped}$", LeftRight, LeftTicks);
