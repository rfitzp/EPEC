import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Stage2/qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] q  = A[1];
real[] x  = A[8];

file    in1 = input("../../Stage2/q95.txt").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

real[] q95 = A1[0];
real[] r95 = A1[1];

pen s = white + dotted + 0.5;
pen s1 = blue;
draw(graph(p,q),s,marker(scale(1.5mm)*polygon(3),s1));
s1 = red;
s = white + dotted + 0.5;	
draw(graph(p,x),s,marker(scale(1.5mm)*polygon(4),s1));

xlimits (0.8,1.0,Crop);

s = dotted + black + 1;
xequals (1.,s);
if (q95[0] > 0.)
  {
    xequals (r95[0],s);
    yequals (q95[0],s);
  }

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/r_a$",BottomTop,LeftTicks);
yaxis("$q$",LeftRight,RightTicks);
