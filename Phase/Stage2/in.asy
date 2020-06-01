import graph;
     
size(500,500,IgnoreAspect);

file    in = input("mirnov.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p     = A[0];
real[] chi1  = A[1];
real[] chi2  = A[2];
real[] chi3  = A[3];
real[] chi4  = A[4];

pen s = blue + 1.5;	
draw(graph(p,chi1),s);
s = red + 1.5;	
draw(graph(p,chi2),s);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\delta B_{p\,{\rm in}} ({\rm G})$",LeftRight,RightTicks);
