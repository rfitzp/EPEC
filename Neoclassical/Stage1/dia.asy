import graph;
     
size(500,500,IgnoreAspect);

file    in = input("diamagnetic.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r      = A[0];
real[] wE     = A[1];
real[] waste  = A[2];
real[] wasti  = A[3];
real[] wastI  = A[4];

pen s = black + dotted + 1.5;	
draw(graph(r,wE),s,marker(scale(1.5mm)*polygon(4)));
s = red + dotted + 1.5;	
draw(graph(r,waste),s,marker(scale(1.5mm)*polygon(4)));
s = blue + dotted + 1.5;	
draw(graph(r,wasti),s,marker(scale(1.5mm)*polygon(4)));
s = green + dotted + 1.5;	
draw(graph(r,wastI),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.4,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$\omega ({\rm krad/s})$",LeftRight,RightTicks);
