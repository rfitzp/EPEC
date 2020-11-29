import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] q   = A[2];

int n = psi.length;

real[] psi1;
real[] q1;

for (int i = 0; i < n; ++i)
{
if (psi[i] > 0.8)
{
psi1.push (psi[i]);
q1.push (q[i]);
}
}

pen s = black + solid + 1.5;	
draw(graph(psi1,q1),s);

s = dotted + black + 1;

xlimits (0.8,1.0,Crop);

yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$n_e(10^{19}\,{\rm m}^{-3})$",LeftRight,RightTicks);
