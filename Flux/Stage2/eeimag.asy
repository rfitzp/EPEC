import graph;

size(500,500,IgnoreAspect);

file    inx = input("../Stage3/E_vector.txt").line();
real[][] Ax = inx.dimension (0,0);
Ax          = transpose(Ax);

real[] EE1   = Ax[2];
real[] EE2   = Ax[4];

file    iny = input("../Stage3/rational.txt").line();
real[][] Ay = iny.dimension (0,0);
Ay          = transpose(Ay);

real[] m = Ay[0];

pen s = black + 1;

draw(graph(m,EE1),s,marker(scale(1.5mm)*polygon(4)));

s = black + dashed;
draw(graph(m,EE2),s,marker(scale(1.5mm)*polygon(3)));

s = black + dotted + 1.5;
yequals (0., s);

//limits ((0.,0.), (2.,2.), Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$m$",BottomTop,LeftTicks);
yaxis("${\rm Im}(\hat{e}_k)$",LeftRight,RightTicks);
