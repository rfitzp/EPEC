import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("diffusivity.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi    = A[0];
real[] chin   = A[1];
real[] chin1  = A[2];
real[] chie   = A[3];
real[] chie1  = A[4];
real[] chip   = A[5];
real[] chip1  = A[6];

real[] chin_1 = chin - 0.5*chin1;
real[] chin_2 = chin + 0.5*chin1;
real[] chie_1 = chie - 0.5*chie1;
real[] chie_2 = chie + 0.5*chie1;
real[] chip_1 = chip - 0.5*chip1;
real[] chip_2 = chip + 0.5*chip1;

int N = psi.length;
real[] x1;
real[] x2;
real[] x3;
for (int j = 0; j < N; ++j)
{
x1.push ( chin_1[j] * fabs(cos(j*pi/2.)) + chin_2[j] * fabs(cos((j+1.)*pi/2.)) );
x2.push ( chie_1[j] * fabs(cos(j*pi/2.)) + chie_2[j] * fabs(cos((j+1.)*pi/2.)) );
x3.push ( chip_1[j] * fabs(cos(j*pi/2.)) + chip_2[j] * fabs(cos((j+1.)*pi/2.)) );
}

pen s = red + solid + 4.0;	
draw(graph(psi,x3),s);
s = green + solid + 4.0;	
draw(graph(psi,x2),s);
s = blue + solid + 4.0;	
draw(graph(psi,x1),s);

s = dotted + black + 1;

limits ((0.,0.), (1.,13.) ,Crop);

//yequals (0., s);

s = dotted + black + 2;

xequals (0.925, s);

pen qq = fontsize(30.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\chi_\perp ({\rm m^{\,2}/s})$",LeftRight,RightTicks);
