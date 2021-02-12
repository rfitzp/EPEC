import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/factor.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    inx = input ("../../../Flux/Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp  = Ax[0];
real psilim = ppp[0];
real psiped = ppp[1];

real[] psi = A[0];
real[] y1  = A[1];
real[] y2  = A[2];
real[] y3  = A[3];
real[] y4  = A[4];

pen s = black + solid + 1.5;	
draw(graph(psi,y1),s,marker(scale(2.0mm)*polygon(4)));
s = red + solid + 1.5;	
draw(graph(psi,y2),s,marker(scale(2.0mm)*polygon(4)));
s = green + solid + 1.5;	
draw(graph(psi,y3),s,marker(scale(2.0mm)*polygon(4)));
s = blue + solid + 1.5;	
draw(graph(psi,y4),s,marker(scale(2.0mm)*polygon(4)));

xlimits (0.85,1.0,Crop);
//ylimits (-40, 60.,Crop);

s = dotted + black + 1.5;
yequals (0.,s);
xequals (psilim, s);
xequals (psiped, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("Factors(1-4) / $10^{19}\,{\rm keV}$",LeftRight,RightTicks);
