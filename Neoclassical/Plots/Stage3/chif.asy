import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/chip.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    inx = input ("../../../Flux/Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp  = Ax[0];
real psilim = ppp[0];
real psiped = ppp[1];

real[] psi    = A[0];
real[] chin   = A[4]+0.01;
real[] chie   = A[3];
real[] chip   = A[2]-0.01;

pen s = red + solid + 4.0;	
draw(graph(psi,chip),s);
draw ((0.7,2.2)--(0.8,2.2),s);

s = blue + solid + 4.0;	
draw(graph(psi,chin),s);
draw ((0.7,1.8)--(0.8,1.8),s);

s = green + solid + 4.0;	
draw(graph(psi,chie),s);
draw ((0.7,2.)--(0.8,2.),s);


s = dotted + black + 1;

limits ((0.,0.), (1.,2.5) ,Crop);

s = dotted + black + 1.5;

xequals (psilim, s);
xequals (psiped, s);

pen pp = fontsize(30.);
label(Label("$\chi_\phi$",(0.83,2.2)));
label(Label("$\chi_E$",(0.83,2.)));
label(Label("$D_\perp$",(0.83,1.8)));

pen qq = fontsize(50.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\chi_\phi,\chi_E,D_\perp ({\rm m^{\,2}/s})$",LeftRight,RightTicks);
