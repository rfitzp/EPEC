import graph;

size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/opt.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

real[] q  = A[1];
real[] pu = A[3];

pen s = black + 1.5;
draw (graph (q, pu), s);

ylimits (0., 2., Crop);

s = black + dotted + 1.5;

yequals (0.5, s);
yequals (1.0, s);
yequals (1.5, s);



pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$q_{95}$",                   BottomTop, RightTicks(Step=0.5,n=5));
yaxis ("${\mit\Delta}/\pi$", LeftRight, LeftTicks(Step = 0.5, n=5));
