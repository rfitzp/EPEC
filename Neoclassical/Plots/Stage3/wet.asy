import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/omega.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[4];
real[] we  = A[7];
real[] wt  = A[8];
real[] wa  = -A[9];
real[] wp  = -A[10];

pen s = black + solid + 1.5;	
draw(graph(psi,we),s,marker(scale(2.0mm)*polygon(4)));
s = red + solid + 1.5;	
draw(graph(psi,wt),s,marker(scale(2.0mm)*polygon(4)));
s = green + solid + 1.5;	
draw(graph(psi,wa),s,marker(scale(2.0mm)*polygon(4)));
s = blue + solid + 1.5;	
draw(graph(psi,wp),s,marker(scale(2.0mm)*polygon(4)));

xlimits (0.85,1.0,Crop);
ylimits (-40, 60.,Crop);

s = dotted + black + 1;
yequals (0.,s);

s = dotted + black + 2;

xequals (0.945, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\omega_E ({\rm krad/s})$",LeftRight,RightTicks);
