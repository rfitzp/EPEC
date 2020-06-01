import graph;

size(500,500,IgnoreAspect);

file    inz = input("../Stage3/parameters.txt").line();
real[][] Az = inz.dimension (0,0);
Az          = transpose(Az);

real[] Nres = Az[9];
int nres   = (int) (Nres[0]);

file    inx = input("../Stage3/E_matrix.txt").line();
real[][] Ax = inx.dimension (0,0);
Ax          = transpose(Ax);

real[] EE   = Ax[2];

file    iny = input("../Stage3/rational.txt").line();
real[][] Ay = iny.dimension (0,0);
Ay          = transpose(Ay);

real[] m = Ay[0];

real[] EEk;
for (int i = 0; i < nres*nres; i += nres+1)
  {
    EEk.push (EE[i]);
  }

pen s = black + 1;

draw(graph(m,EEk),s,marker(scale(1.5mm)*polygon(4)));

s = dotted + 1.5;
draw ((0.,0.)--(m[m.length-1],-2.*m[m.length-1]), s);

s = white + dotted + 0.;
yequals (0., s);
xequals (0., s);

//limits ((0.,0.), (2.,2.), Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$m$",BottomTop,LeftTicks);
yaxis("$E_{kk}$",LeftRight,RightTicks);
