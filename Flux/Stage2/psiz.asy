import graph;
import palette;
import contour;
   
file    in = input ("PsiZ.txt").line();
real[][] a = in.dimension (0, 0);
  
file    in2 = input("../Stage1/Box.txt").line();
real[][] A2 = in2.dimension (0, 0);
A2          = transpose (A2);
real[] aa   = A2[0];
real rmin   = (real) aa[0];
real[] bb   = A2[1];
real zmin   = (real) bb[0]; 
real[] cc   = A2[2];
real rmax   = (real) cc[0]; 
real[] dd   = A2[3];
real zmax   = (real) dd[0]; 
real zsize  = 500;
real rsize  = 500*(rmax-rmin)/(zmax-zmin);

file    in4 = input("../Stage1/Boundary.txt").line();
real[][] A4 = in4.dimension (0, 0);
A4          = transpose (A4);
real[] rb   = A4[0];
real[] zb   = A4[1];

file    in3 = input("../Stage1/Axis.txt").line();
real[][] A3 = in3.dimension (0, 0);
A3          = transpose (A3);
real[] xx   = A3[0];
real rw     = (real) xx[0];
real[] yy   = A3[1];
real zw     = (real) yy[0]; 

size (rsize, zsize, (rmin, zmin), (rmax, zmax));

pen[] Palette = Rainbow();
bounds range  = image (a, Full, (rmin, zmin), (rmax, zmax), Palette);

pen s = white + 1.5;
draw (graph (rb, zb), s);

real[] cvals = uniform (range.min, range.max, 40);
draw (contour (a, (rmin, zmin), (rmax, zmax), cvals));
pen s = red + 1.5;
draw (contour (a, (rmin, zmin), (rmax, zmax), new real[] {0.}),s);

s = white;
dot ((rw, zw), s);

pen q = fontsize (20.);
defaultpen (q);
xaxis ("$R/R_0$",  Bottom, RightTicks, above=true);
xaxis (Top, NoTicks, above=true);
yaxis ("$Z/Z_0$", Right, RightTicks, above=true);
yaxis (Left, NoTicks, above=true);
