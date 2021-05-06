import graph;

size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/opt.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

real[] t  = A[0];
real[] pu = A[5];
real[] pm = A[6];
real[] pl = A[7];

pen s = red + 1.5;
draw (graph (t, pu), s);

s = green + 1.5;
draw (graph (t, pm), s);

s = blue + 1.5;
draw (graph (t, pl), s);

ylimits (-2., 2., Crop);

s = black + dotted + 1.5;

yequals (-0.5, s);
yequals (-1.0, s);
yequals (-1.5, s);
yequals (0., s);
yequals (0.5, s);
yequals (1.0, s);
yequals (1.5, s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",                   BottomTop, RightTicks);
yaxis ("${\mit\Delta}_{L,U,M}/\pi$", LeftRight, LeftTicks(Step = 0.5, n=5));