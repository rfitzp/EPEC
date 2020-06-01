import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Flux/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi    = A[0];
real[] r      = A[1];
real[] psir   = A[4];

pen s = black + dotted + 1.5;	
draw(graph(r,psi),s,marker(scale(1.5mm)*polygon(4)));

xlimits (0.,1.0,Crop);

s = dotted + black + 1;
//yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$\hat{r}$",BottomTop,LeftTicks);
yaxis("${\mit\Psi}_N$",LeftRight,RightTicks);
