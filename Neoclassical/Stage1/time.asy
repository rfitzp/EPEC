import graph;
     
size(500,500,IgnoreAspect);

file    in = input("timescale.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r      = A[0];
real[] ta     = A[1];
real[] tr     = A[2];
real[] tm     = A[3];
real[] tth    = A[4];

pen s = black + dotted + 1.5;	
draw(graph(r,ta),s,marker(scale(1.5mm)*polygon(4)));
s = red + dotted + 1.5;	
draw(graph(r,tr),s,marker(scale(1.5mm)*polygon(4)));
s = blue + dotted + 1.5;	
draw(graph(r,tm),s,marker(scale(1.5mm)*polygon(4)));
s = green + dotted + 1.5;	
draw(graph(r,tth),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$\log_{10}[\tau({\rm s})]$",LeftRight,RightTicks);
