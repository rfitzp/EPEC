import graph;
     
size(500,500,IgnoreAspect);

file    in = input("qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] q  = A[1];
real[] x  = A[5];

pen s = solid + 0.;
yequals (0.);
//s = red + dotted + 0.2;	
//draw(graph(p,q),s,marker(scale(1.5mm)*polygon(4)));
s = green + dotted + 0.2;	
draw(graph(p,x),s,MarkFill[0]);

s = dotted + black + 1;
xequals (1.,s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/r_a$",BottomTop,LeftTicks);
yaxis("$p$",LeftRight,RightTicks);
