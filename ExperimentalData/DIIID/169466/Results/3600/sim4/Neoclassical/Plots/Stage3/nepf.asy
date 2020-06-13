import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi  = A[0];
real[] r    = A[1];
real[] q    = A[3];
real[] dpdr = A[19];
real[] qq   = q/dpdr;

pen s = black + solid + 1.5;	
draw(graph(r,qq),s);

s = dotted + black + 1;

xlimits (0.,1.0,Crop);

yequals (0., s);

pen qqx = fontsize(25.);
defaultpen (qqx);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$dn_e/d{\rm \Psi}_N(10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
import graph;
     
