import graph;

file    inz = input("nres.txt").line();
real[][] Az = inz.dimension (0,0);
Az          = transpose(Az);

real[] Nres = Az[0];
real nres   = Nres[0];

file    inx = input("Rst.txt").line();
real[][] Ax = inx.dimension (0,0);
Ax          = transpose(Ax);

file    iny = input("Zst.txt").line();
real[][] Ay = iny.dimension (0,0);
Ay          = transpose(Ay);
  
file    in2 = input("../Stage1/Box.txt").line();
real[][] A2 = in2.dimension (0,0);
A2          = transpose(A2);
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
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] rb   = A4[0];
real[] zb   = A4[1];

file    in3 = input("../Stage1/Axis.txt").line();
real[][] A3 = in3.dimension (0,0);
A3          = transpose(A3);
real[] xx   = A3[0];
real rw     = (real) xx[0];
real[] yy   = A3[1];
real zw     = (real) yy[0]; 

size (rsize, zsize, (rmin, zmin), (rmax, zmax));

pen s = black+3.;
draw ((rmin,zmin)--(rmax,zmin)--(rmax,zmax)--(rmin,zmax)--(rmin,zmin),s);

s = blue+0.25;
draw(graph(rb,zb),s);

s = blue;
dot ((rw, zw), s);

s = black + 0.25;
for (int i = 0; i < nres; ++i)
  {
    real[] Rst = Ax[i+1];
    real[] Zst = Ay[i+1];

    //draw(graph(Rst,Zst),s,marker(scale(0.5mm)*polygon(3)));//
    draw(graph(Rst,Zst),s);
  }

pen q = fontsize(20.);
defaultpen (q);
xaxis("$R/R_0$",  Bottom, RightTicks, above=true);
xaxis(Top, NoTicks, above=true);
yaxis("$Z/R_0$", Right, RightTicks, above=true);
yaxis(Left, NoTicks, above=true);
