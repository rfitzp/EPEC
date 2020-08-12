import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("profiles.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] ne1 = A[1];
real[] ne2 = A[2];
real[] ne3 = A[3];
real[] te1 = A[4];
real[] te2 = A[5];
real[] te3 = A[6];
real[] ti1 = A[7];
real[] ti2 = A[8];
real[] ti3 = A[9];
real[] ni1 = A[10]*10.;
real[] ni2 = A[11]*10.;
real[] ni3 = A[12]*10.;

int N = psi.length;
real[] x1;
real[] x2;
real[] x3;
real[] x4;
for (int j = 0; j < N; ++j)
{
x1.push ( ne1[j] * fabs(cos(j*pi/2.)) + ne3[j] * fabs(cos((j+1.)*pi/2.)) );
x2.push ( te1[j] * fabs(cos(j*pi/2.)) + te3[j] * fabs(cos((j+1.)*pi/2.)) );
x3.push ( ti1[j] * fabs(cos(j*pi/2.)) + ti3[j] * fabs(cos((j+1.)*pi/2.)) );
x4.push ( ni1[j] * fabs(cos(j*pi/2.)) + ni3[j] * fabs(cos((j+1.)*pi/2.)) );
}

pen s = red + solid + 4.;	
draw(graph(psi,x1),s);


s = green + solid + 4.0;	
draw(graph(psi,x2),s);

s = blue + solid + 4.0;	
draw(graph(psi,x3),s);

s = cyan + solid + 4.0;	
draw(graph(psi,x4),s);

s = dotted + black + 1;

limits ((0.,0.),(1.,6.5),Crop);

//yequals (0., s);

s = dotted + black + 2;

xequals (0.925, s);

pen qq = fontsize(30.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$n_{e,I}$, $T_{e,i}$",LeftRight,RightTicks);
