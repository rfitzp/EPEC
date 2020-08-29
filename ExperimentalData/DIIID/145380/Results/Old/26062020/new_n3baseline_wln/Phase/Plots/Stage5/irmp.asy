import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/rmp.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] t  = A[0];
real[] i  = A[1];
real[] p  = A[2];

pen s = black + 1.5;	
draw(graph(t,i),s);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$I_{\rm rmp} ({\rm kA})$",LeftRight,RightTicks);