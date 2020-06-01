import graph;
     
size(500,500,IgnoreAspect);

file    in = input("qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] q  = A[1];
real[] x  = A[8];

file    in1 = input("q95.txt").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

real[] q95 = A1[0];
real[] r95 = A1[1];

pen s = solid + 0.;
yequals (0.);

s      = red + dotted + 0.2;
pen s1 = blue;
draw(graph(p,q),s,marker(scale(0.5mm)*polygon(3), s1));
s  = green + dotted + 0.2;
s1 = red;
draw(graph(p,x),s,marker(scale(0.5mm)*polygon(3), s1));

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
