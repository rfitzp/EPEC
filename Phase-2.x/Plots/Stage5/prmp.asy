import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/rmp.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] t  = A[0];
real[] i  = A[1];
real[] p  = A[2];

pen s = solid + red + 4.;	
draw(graph(t,p),s);

s = dotted + black + 1.0;
yequals (0., s);

ylimits (0., 2., Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("${\mit\Delta}\phi_{\rm rmp}/\pi$",LeftRight,RightTicks);
