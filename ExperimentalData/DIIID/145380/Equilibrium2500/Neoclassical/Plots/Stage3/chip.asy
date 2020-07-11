import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/chip.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] q   = A[2];
real[] x   = A[3];
real[] y   = A[4];

pen s = red + solid + 3.0;	
draw(graph(psi,q),s);
s = green + solid + 3.0;	
draw(graph(psi,x),s);
s = blue + solid + 3.0;	
draw(graph(psi,y),s);

s = dotted + black + 1;

limits ((0.,0.), (1.,13.) ,Crop);

//yequals (0., s);

s = dotted + black + 2;

xequals (0.925, s);

pen qq = fontsize(30.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\chi_\perp ({\rm m^{\,2}/s})$",LeftRight,RightTicks);
