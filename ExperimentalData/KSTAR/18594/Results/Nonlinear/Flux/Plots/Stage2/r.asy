import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage2/qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] x  = A[9];
real[] y  = A[10];

pen s = green + 1.;
draw(graph(p,x),s);
s = blue + 1.;
draw(graph(p,y),s);

xlimits (0., 1., Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/r_a$",BottomTop,LeftTicks);
yaxis("$R/R_0$",LeftRight,RightTicks);
