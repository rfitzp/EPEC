import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    inx = input ("../../../Flux/Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp  = Ax[0];
real psilim = ppp[2];
real psiped = ppp[1];
     
real[] psi = A[0];
real[] ne1 = A[28];

pen s  = red + dotted + 0.2;
pen s1 = blue;
draw(graph(psi,ne1),s, marker(scale(0.5mm)*polygon(3), s1));

xlimits (0.85,1.0,Crop);

s = dotted + 1.5 + black;
yequals (0., s);
xequals (psilim, s);
xequals (psiped, s);

pen qqx = fontsize(25.);
defaultpen (qqx);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$d^{\,2}n_i/d{\rm \Psi}_N^{\,2}(10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
import graph;
     
