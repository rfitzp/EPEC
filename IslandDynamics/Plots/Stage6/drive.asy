import graph;

size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/opt.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

real[] t  = A[0];
real[] pu = A[8];
real[] pm = A[9];

pen s = red + 1.5;
draw (graph (t, pm), s);

s = blue + 1.5;
draw (graph (t, pu), s);

s = black + 1 + dotted;

yequals (0., s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",                   BottomTop, RightTicks);
//yaxis ("$\chi_k$", LeftRight, LeftTicks(Step=5e-5,n=5));
yaxis ("$\chi_k$", LeftRight, LeftTicks);