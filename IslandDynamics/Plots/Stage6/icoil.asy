import graph;

size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/opt.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

real[] t  = A[0];
real[] iu = A[2];
real[] im = A[3];
real[] il = A[4];

pen s = red + 1.;
draw (graph (t, iu), s);

s = green + 1.;
draw (graph (t, im), s);

s = blue + 1.;
draw (graph (t, il), s);

s = black + 1 + dotted;

yequals (0., s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",                BottomTop, RightTicks);
yaxis ("$I_{L,U,M}({\rm kA/turn})$", LeftRight, LeftTicks);