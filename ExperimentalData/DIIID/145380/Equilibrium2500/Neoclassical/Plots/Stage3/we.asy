import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("profiles.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] w1  = A[13];
real[] w2  = A[14];
real[] w3  = A[15];

int N = psi.length;
real[] x1;
for (int j = 0; j < N; ++j)
{
x1.push ( w1[j] * fabs(cos(j*pi/2.)) + w3[j] * fabs(cos((j+1.)*pi/2.)) );
}

pen s = red + solid + 4.0;	
draw(graph(psi,x1),s);

s = dotted + black + 1.0;
yequals (0.,s);

s = dotted + black + 2;

xequals (0.925, s);

xlimits (0.,1.0,Crop);

pen qq = fontsize(30.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\omega_E ({\rm krad/s})$",LeftRight,RightTicks);
