import graph;
     
size(500,500,IgnoreAspect);

file    in = input("profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] p   = A[4];
real[] pe  = A[5];
real[] pi  = A[6];

pen s;

s = black + 1.;	
draw(graph(r,p),s);
s = red + 1.;	
draw(graph(r,pe),s);
s = blue + 1.;	
draw(graph(r,pi),s);

s = dotted + black + 1;
//yequals (0.,s);

xlimits (0.,1.0,Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$p$",LeftRight,RightTicks);
