import graph;
     
size(750,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/bootstrap.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r   = A[0];
real[] p   = A[1];
real[] abe = A[2];
real[] abi = A[3];
real[] ac  = A[4];
real[] ap  = A[5];

real[] abc = abe + abi + ac;

ap = ap*1.e4;

pen s = red + dotted + 1.5;	
draw(graph(p,abe),s,marker(scale(1.5mm)*polygon(4)));
s = blue + dotted + 1.5;	
draw(graph(p,abi),s,marker(scale(1.5mm)*polygon(4)));
s = green + dotted + 1.5;	
draw(graph(p,ac),s,marker(scale(1.5mm)*polygon(4)));
s = black + dotted + 1.5;	
draw(graph(p,abc),s,marker(scale(1.5mm)*polygon(4)));
s = cyan + dotted + 1.5;	
draw(graph(p,ap),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.2,1.0,Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\alpha_b$, $\alpha_c$, $10^4\,\alpha_p$",LeftRight,RightTicks);
