import graph;
import palette;
import contour;
     
file    in2 = input ("../../Outputs/Stage1/Box.txt").line();
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
real rsize  = 500 * (rmax - rmin) /(zmax - zmin);

size (rsize, zsize, (rmin, zmin), (rmax, zmax));

file    in1 = input ("../../Outputs/Stage1/Boundary.txt").line();
real[][] A1 = in1.dimension (0, 0);
A1          = transpose (A1);
     
real[] r  = A1[0];
real[] z  = A1[1];

file    in3 = input ("../../Outputs/Stage1/Limiter.txt").line();
real[][] A3 = in3.dimension (0, 0);
A3          = transpose (A3);
     
real[] rl  = A3[0];
real[] zl  = A3[1];

file    in4 = input ("../../Outputs/Stage1/Axis.txt").line();
real[][] A4 = in4.dimension (0, 0);
A4          = transpose (A4);
     
real[] Raxis = A4[0];
real[] Zaxis = A4[1];
real raxis   = (real) Raxis[0];
real zaxis   = (real) Zaxis[0];

file    in = input ("../../Outputs/Stage1/Psi.txt").line();
real[][] a = in.dimension (0, 0);
     
pen[] Palette = Rainbow ();
bounds range  = image (a, Full, (rmin, zmin), (rmax, zmax), Palette);

real[] cvals = uniform (0., range.max, 20);
draw (contour (a, (rmin, zmin), (rmax, zmax), cvals));

cvals = uniform (range.min, 0., 40);
draw (contour (a, (rmin, zmin), (rmax, zmax), cvals));

pen s = solid + white + 2.5;
dot ((raxis, zaxis), s);

s = solid + white + 1.5;
draw (graph (r, z), s);

s = solid + black + 2.5;

draw (graph (rl, zl), s);

pen q = fontsize (20.);
defaultpen (q);
//xaxis ("$R/R_0$", Bottom, RightTicks(new real[] {0.75,1.,1.25, 1.5}), above = true);
xaxis ("$R/R_0$", Bottom, RightTicks, above = true);
xaxis (Top,  NoTicks, above = true);
yaxis ("$Z/R_0$", Right,  RightTicks, above = true);
yaxis (Left, NoTicks, above = true);



