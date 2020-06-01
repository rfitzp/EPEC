import graph;
     
size(500,500,IgnoreAspect);

file    in = input("phi.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p     = A[0];

int n = getint ("n ?");

pen s = black + 1.5;	
draw(graph(p,A[n]),s);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\varphi/\pi$",LeftRight,RightTicks);
