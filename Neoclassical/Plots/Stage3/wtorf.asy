import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] q   = A[20];

pen s = red + solid + 4.;	
draw(graph(psi,q),s);

s = dotted + black + 1;
yequals (0.,s);

s = dotted + black + 2;

xequals (0.945, s);

xlimits (0.,1.0,Crop);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\omega_{\phi\,I} ({\rm krad/s})$",LeftRight,RightTicks);
