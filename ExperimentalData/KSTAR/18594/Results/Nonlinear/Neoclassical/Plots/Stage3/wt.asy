import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/omega.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r    = A[0];
real[] web  = A[6];
real[] webI = A[7];
real[] wtor = A[8];
real[] wast = A[9];
real[] wthe = A[10];

pen s = black + dotted + 1.5;	
draw(graph(r,wtor),s,marker(scale(1.5mm)*polygon(4)));
s = red + dotted + 1.5;	
draw(graph(r,webI),s,marker(scale(1.5mm)*polygon(4)));
s = green + dotted + 1.5;	
draw(graph(r,wast),s,marker(scale(1.5mm)*polygon(4)));
s = blue + dotted + 1.5;	
draw(graph(r,wthe),s,marker(scale(1.5mm)*polygon(4)));

s = yellow + dotted + 1.5;	
draw(graph(r,web),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.4,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$\omega_\phi ({\rm krad/s})$",LeftRight,RightTicks);
