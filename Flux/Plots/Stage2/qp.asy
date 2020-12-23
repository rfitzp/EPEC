import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage2/qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[3];
real[] q  = A[1];

real[] psi = 1. - p;

pen s  = red + dotted + 0.2;
pen s1 = blue;
draw(graph(psi,q),s,marker(scale(0.5mm)*polygon(3), s1));

xlimits (0.85, 0.99, Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$q$",LeftRight,RightTicks);
