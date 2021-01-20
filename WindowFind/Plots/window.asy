import graph;
import palette;
import contour;

size (1000, 500, IgnoreAspect);

file    in  = input ("../Outputs/window.txt").line();
real[][] A  = in.dimension (0, 0);
A           = transpose (A);
file    in1 = input ("../Outputs/limits.txt").line();
real[][] A1 = in1.dimension (0, 0);
A1          = transpose (A1);

real xmin = A1[0][0];
real xmax = A1[1][0];
real ymin = A1[2][0];
real ymax = A1[3][0];

pen[] Palette = BWRainbow ();
bounds range  = image (A, Full, (xmin, ymin), (xmax, ymax), Palette);

pen s = 1.5 + black;
real[] cvals = {0.05};
draw (contour (A, (xmin, ymin), (xmax, ymax), cvals), s);

pen q = fontsize (20.);
defaultpen (q);
xaxis ("$q_{95}$", Bottom, RightTicks (n=5), above=false);
xaxis (Top,  NoTicks);
yaxis ("$I_{\rm rmp}({\rm kA})$", Left, LeftTicks);
yaxis (Right, NoTicks);

palette ("$-{\mit\Delta}P/P$", range, point(NW) + (0, 0.02), point(NE) + (0, 0.1), Top, Palette);

