import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/omega.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    inx = input ("../../../Flux/Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp  = Ax[0];
real psilim = ppp[2];
real psiped = ppp[1];

real[] psi = A[4];
real[] we  = A[6];
real[] wt  = A[7];
real[] wa  = -A[8];
real[] wp  = -A[9];

pen s = black + solid + 1.5;	
draw(graph(psi,we),s,marker(scale(2.0mm)*polygon(4)));
s = red + dashed + 1.5;	
draw(graph(psi,wt),s,marker(scale(2.0mm)*polygon(4)));
s = green + solid + 1.5;	
draw(graph(psi,wa),s,marker(scale(2.0mm)*polygon(4)));
s = blue + solid + 1.5;	
draw(graph(psi,wp),s,marker(scale(2.0mm)*polygon(4)));

xlimits (0.85,1.0,Crop);
//ylimits (-40, 60.,Crop);

s = dotted + black + 1;
yequals (0.,s);

s = dotted + black + 1.5;
xequals (psilim, s);
xequals (psiped, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\omega_E ({\rm krad/s})$",LeftRight,RightTicks);
