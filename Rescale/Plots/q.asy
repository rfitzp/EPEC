import graph;
     
size (500, 500, IgnoreAspect);

file    in = input ("../Outputs/Profiles.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] p  = A[0];
real[] g  = A[5];

file    in1 = input ("../Outputs/Profiles1.txt").line();
real[][] A1 = in1.dimension (0, 0);
A1          = transpose (A1);
     
real[] p1  = A1[0];
real[] g1  = A1[5];

pen s = solid + 1.5;
draw (graph (p, g), s, marker(scale(1.5mm)*polygon(4)));
s = solid + blue + 1.5;
pen ss = solid + blue + 0.5;
draw (graph (p1, g1), s, marker(scale(1.5mm)*polygon(4),ss));

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("$\Psi_N$", BottomTop, LeftTicks);
yaxis("$q$",      LeftRight, RightTicks);
