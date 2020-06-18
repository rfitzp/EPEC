import graph;
     
size (750, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/q.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] q0 = A[0];
real[] q  = A[1];
real[] qa = A[2];
real[] t  = A[4];

pen s = black + dotted + 1.;

draw (graph (t, q),  s, marker (scale (1.mm)*polygon (10), green));

s = dotted + black + 1;

//ylimits (0.,5., Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",           BottomTop, LeftTicks);
yaxis ("$q_{95}$", LeftRight, RightTicks);
