import graph;
     
size(500,500,IgnoreAspect);

file    in = input("profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] q   = A[18];

pen s = black + solid + 1.5;	
draw(graph(r,q),s);

s = dotted + black + 1;

xlimits (0.,1.0,Crop);

//yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$n_b(10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
