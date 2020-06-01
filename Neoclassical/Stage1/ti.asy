import graph;
     
size(500,500,IgnoreAspect);

file    in = input("profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] q   = A[8];

pen s = black + 1.5;
draw(graph(r,q),s);

s = dotted + black + 1;
//yequals (0.,s);

xlimits (0.,1.0,Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$T_i ({\rm keV})$",LeftRight,RightTicks);
