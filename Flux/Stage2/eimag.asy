import graph;
import palette;
import contour;

size (500, 500, IgnoreAspect);


file    inz = input("../Stage3/parameters.txt").line();
real[][] Az = inz.dimension (0,0);
Az          = transpose(Az);

real[] Nres = Az[9];
int nres   = (int) (Nres[0]);

file    inx = input("../Stage3/E_matrix.txt").line();
real[][] Ax = inx.dimension (0,0);
Ax          = transpose(Ax);

real[] m1   = Ax[0];
real[] m2   = Ax[1];
real[] EE   = Ax[3];

pair[] z;
real[] f;

for (int i = 0; i < nres; i += 1)
  {
    for (int j = 0; j < nres; j += 1)
      {
	z.push ((m1[i+nres*j], m2[i+nres*j]));
	f.push (fabs(EE[i+nres*j]));
      }
  }

pen[] Palette = BWRainbow ();
bounds range  = image (z, f, Full, Palette);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("$m$", BottomTop, LeftTicks);
yaxis("$m$", LeftRight, RightTicks);

qq = fontsize (15.);
defaultpen (qq);
pen Tickpen = black + 1;
pen tickpen = gray  + 0.5*linewidth(currentpen);
palette (range, point(NW) + (0, 0.5), point(NE) + (0, 1), Top, Palette,
        PaletteTicks ("$%#.2f$", N = 10,n = 0, Tickpen, tickpen));
