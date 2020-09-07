import graph;
import palette;
import contour;

size (500,500,IgnoreAspect);

real gamma = getreal ("scale ");
real skip  = getint  ("skip");
real t1    = getreal ("T_start ");
real t2    = getreal ("T_end ");
write (" ");

file    in = input ("../../Outputs/Stage5/mirnov.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

real[] t = A[0];
real[] s = A[3];
real[] c = A[4];

file    in1 = input ("../../Outputs/Stage5/ntor.txt").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

real[] n = A1[0];

real ntor = n[0];

int l  = t.length;

pair[] z;
real[] f;
for (int i = 0; i < l; i += 1)
  {
    for (int j = 0; j < 360; j += 1)
      {
	real p = j * pi /180.;
	real b = c[i]*cos(ntor*p) + s[i]*sin(ntor*p);

	if (t[i] > t1 && t[i] < t2 && i%skip == 0)
	  {
	    z.push ((t[i], p/pi));
	    if (b > 0.)
	      f.push (b**gamma);
	    else
	      f.push (-(-b)**gamma);
	  }
      }
  }

pen[] Palette = BWRainbow();
//pen[] Palette = Grayscale();
bounds range = image (z, f, Full, Palette);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\phi/\pi$",LeftRight,RightTicks);

qq = fontsize (15.);
defaultpen (qq);
pen Tickpen = black + 1;
pen tickpen = gray  + 0.5*linewidth(currentpen);
palette (range, point(NW)+(0,0.2), point(NE)+(0,0.45), Top, Palette,
        PaletteTicks ("$%#.2f$", N=10,n=0,Tickpen,tickpen));
