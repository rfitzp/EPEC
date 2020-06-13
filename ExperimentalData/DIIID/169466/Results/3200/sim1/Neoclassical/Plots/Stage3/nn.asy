import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/diamagnetic.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r      = A[0];
real[] wE     = A[1];
real[] waste  = A[2];
real[] wasti  = A[3];
real[] wastI  = A[6];

pen s = black + dotted + 1.5;	
s = black + dotted + 1.5;	
draw(graph(r,wastI),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.4,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$n_n(10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
