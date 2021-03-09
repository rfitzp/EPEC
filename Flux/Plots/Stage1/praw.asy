import graph;
     
size (500, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage1/RawProfiles.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

file    inx = input ("../../Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp = Ax[0];
real psilim = ppp[0];
real psiped = ppp[1];
     
real[] p  = A[0];
real[] g  = A[2]/1.e4;

pen s = solid + 1.5;
draw (graph (p, g), s, marker(scale(1.5mm)*polygon(4)));

s = dotted + 1.5;
xequals (psilim, s);
xequals (psiped, s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$", BottomTop, LeftTicks);
yaxis("$P (10^{\,4}\,{\rm Pa})$",      LeftRight, RightTicks);
