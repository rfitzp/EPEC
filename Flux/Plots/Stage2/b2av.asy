import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage2/qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] x  = A[15];

pen s = black+ 1.;
draw(graph(p,x),s);

xlimits (0., 1., Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/r_a$",BottomTop,LeftTicks);
yaxis("$\langle B^{\,2}\rangle/B_0^{\,2}$",LeftRight,RightTicks);
