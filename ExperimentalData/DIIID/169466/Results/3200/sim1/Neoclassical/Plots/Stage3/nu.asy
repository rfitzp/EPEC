import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/neoclassical.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r      = A[0];
real[] nPe    = A[2];
real[] nPi    = A[3];
real[] nPI    = A[4];
real[] nPSe   = A[5];
real[] nPSi   = A[6];
real[] nPSI   = A[7];

pen s = red + dotted + 1.5;	
draw(graph(r,nPe),s,marker(scale(1.5mm)*polygon(4)));
s = blue + dotted + 1.5;	
draw(graph(r,nPi),s,marker(scale(1.5mm)*polygon(4)));
s = green + dotted + 1.5;	
draw(graph(r,nPI),s,marker(scale(1.5mm)*polygon(4)));

s = red + dotted + 1.5;	
draw(graph(r,nPSe),s,marker(scale(1.5mm)*polygon(3)));
s = blue + dotted + 1.5;	
draw(graph(r,nPSi),s,marker(scale(1.5mm)*polygon(3)));
s = green + dotted + 1.5;	
draw(graph(r,nPSI),s,marker(scale(1.5mm)*polygon(3)));

xlimits (0.,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$\log_{10}(\nu_\ast)$",LeftRight,RightTicks);
