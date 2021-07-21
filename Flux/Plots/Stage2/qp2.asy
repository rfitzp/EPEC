import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage2/qr.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    inx = input ("../../Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp = Ax[0];
real psilim = ppp[0];
real psiped = ppp[1];
         
real[] p  = A[3];
real[] q  = A[16];

for (int i = 0; i < q.length; ++i)
{
if (q[i]>0)
q[i] = q[i]**0.25;
else
q[i] = -(-q[i])**0.25;
}

real[] psi = 1. - p;

pen s  = red + dotted + 0.2;
pen s1 = blue;
draw(graph(psi,q),s,marker(scale(0.5mm)*polygon(3), s1));

xlimits (0.85, 1., Crop);

s = dotted + 1.5 + black;
xequals (psilim, s);
xequals (psiped, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$|d^2q/d{\mit\Psi}_N^2|^{1/4}$",LeftRight,RightTicks);
