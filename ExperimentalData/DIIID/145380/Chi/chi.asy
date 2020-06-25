import graph;
     
size(500,500,IgnoreAspect);

file    in = input("c145380.2500").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] chi = A[1];

pen s = black + solid + 1.5;	
draw(graph(psi,chi), s, marker (scale(2mm)*polygon(3), filltype=Fill, red));

limits ((0.5, 0.), (1., 12.), Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\chi_\perp ({\rm m^{\,2}/s})$",LeftRight,RightTicks);
