import graph;
     
size (500, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage1/RawProfiles.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] p  = A[0];
real[] g  = A[2];

pen s = solid + 1.5;
draw (graph (p, g), s, marker(scale(1.5mm)*polygon(4)));

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("$\psi_N$", BottomTop, LeftTicks);
yaxis("$P$",      LeftRight, RightTicks);
