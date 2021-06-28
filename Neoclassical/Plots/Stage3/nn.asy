import graph;
     
size(750,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] nn  = A[21];

pen s = black + 1.;
draw(graph(r,nn),s);

xlimits (0.4,1.0,Crop);

s = black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$n_n(10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
