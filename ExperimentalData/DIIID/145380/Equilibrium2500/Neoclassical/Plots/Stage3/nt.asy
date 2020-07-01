import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] ne  = A[2];
real[] te  = A[4];
real[] ti  = A[8];
real[] nii = A[10]*10.;

pen s = green + solid + 1.5;	
draw(graph(psi,ne),s);
s = red + solid + 1.5;	
draw(graph(psi,te),s);
s = blue + solid + 1.5;	
draw(graph(psi,ti),s);
s = cyan + solid + 1.5;	
draw(graph(psi,nii),s);

s = dotted + black + 1;

limits ((0.,0.),(1.,6.5),Crop);

//yequals (0., s);

s = dotted + black + 2;

xequals (0.925, s);

pen qq = fontsize(30.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$n_{e,I}$, $T_{e,i}$",LeftRight,RightTicks);
