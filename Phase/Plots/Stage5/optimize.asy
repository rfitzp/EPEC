
import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/optimize.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] t  = A[4];
real[] i  = A[7];
real[] p  = A[8];
real[] q  = i/p;

pen s = solid + red + 3.;	
//draw(graph(t,p),s);
s = solid + blue + 3.;	
draw(graph(t,q),s);

//ylimits (0., 2.5, Crop);

s = dotted + black + 1.;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Delta}/\pi$",BottomTop,LeftTicks(Size=0.1));
yaxis("$\chi$",LeftRight,RightTicks);
