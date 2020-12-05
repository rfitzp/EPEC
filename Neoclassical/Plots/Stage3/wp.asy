import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/omega.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r    = A[4];
real[] web  = A[6];
real[] webI = A[7];
real[] wtor = A[8];
real[] wast = A[9];
real[] wthe = A[10];
real[] wp = wtor - web - wast;

pen s = black + dotted + 1.5;	
draw(graph(r,wp),s,marker(scale(1.5mm)*polygon(4)));
s = red + dotted + 1.5;	
draw(graph(r,wthe),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.4,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\omega_\theta ({\rm krad/s})$",LeftRight,RightTicks);
