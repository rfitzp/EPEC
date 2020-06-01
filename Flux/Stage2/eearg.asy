import graph;

size(500,500,IgnoreAspect);

file    inx = input("../Stage3/E_vector.txt").line();
real[][] Ax = inx.dimension (0,0);
Ax          = transpose(Ax);

real[] EE0   = Ax[1];
real[] EE1   = Ax[2];
real[] EE2   = Ax[3];
real[] EE3   = Ax[4];

real[] x;
real[] y;
for (int i = 0; i < EE0.length; i += 1)
  {
    x.push (atan2(EE1[i],EE0[i])/pi);
    y.push (atan2(EE3[i],EE2[i])/pi);
  }

file    iny = input("../Stage3/rational.txt").line();
real[][] Ay = iny.dimension (0,0);
Ay          = transpose(Ay);

real[] m = Ay[0];

pen s = black + 1;

draw(graph(m,x),s,marker(scale(1.5mm)*polygon(4)));

s = black + dashed;
draw(graph(m,y),s,marker(scale(1.5mm)*polygon(3)));

s = black + dotted + 1.5;
yequals (0., s);

//limits ((0.,0.), (2.,2.), Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$m$",BottomTop,LeftTicks);
yaxis("${\rm arg}(\hat{e}_k)/\pi$",LeftRight,RightTicks);
