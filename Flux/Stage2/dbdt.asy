import graph;

size(500,500,IgnoreAspect);

file    inz = input("nres.out").line();
real[][] Az = inz.dimension (0,0);
Az          = transpose(Az);

real[] Nres = Az[0];
real nres   = Nres[0];

file    inx = input("Cnc.out").line();
real[][] Ax = inx.dimension (0,0);
Ax          = transpose(Ax);

file    iny = input("Znc.out").line();
real[][] Ay = iny.dimension (0,0);
Ay          = transpose(Ay);
  
pen s = black + 1;
real[] theta = Ax[0];
for (int i = 0; i < nres; ++i)
  {
    real[] Rst = Ax[i+1];
    real[] Zst = Ay[i+1];

    //draw(graph(Rst,Zst),s,marker(scale(0.5mm)*polygon(3)));//
    draw(graph(theta,Rst),s);
  }

s = dotted + 1.5;
yequals (0., s);
xequals (1., s);

//limits ((0.,0.), (2.,2.), Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("${\mit\Theta}$",BottomTop,LeftTicks);
yaxis("$dB/d{\mit\Theta}$",LeftRight,RightTicks);
