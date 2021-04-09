import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/optimize.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[4];
real[] x  = A[7];
real[] y  = A[8];

real xm = max(x);
real ym = max(y);
real[] xx = x/xm;
real[] yy = y/ym;

pen s = solid + red + 3.;	
draw (graph (p, xx), s);
s = solid + blue + 3.;	
draw (graph (p, yy), s);

s = dotted + black + 1.;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Delta}/\pi$",BottomTop,LeftTicks);
yaxis("$\chi$",LeftRight,RightTicks);
