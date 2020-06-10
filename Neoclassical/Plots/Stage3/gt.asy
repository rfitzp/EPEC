import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/neoclassical.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r      = A[0];
real[] gt     = A[1];

pen s = black + dotted + 1.5;	
draw(graph(r,gt),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.,1.0,Crop);

s = dotted + black + 1;
//yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$g_t$",LeftRight,RightTicks);