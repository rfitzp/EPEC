import graph;
     
size(500,500,IgnoreAspect);

file    in = input("mirnov.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] t     = A[0];
real[] bcos  = A[1];
real[] bsin  = A[2];

int I = t.length;
int J = 360;

pair[][] z;
for (int i = 0; i < I; ++i)
  {
    for (int j = 0; j < J; ++j)
      {
	z[i][j] = (t[i], real (j) * pi /180.);
      }
  }

pen s = blue + 1.5;	
draw(graph(p,chi1),s);
s = red + 1.5;	
draw(graph(p,chi2),s);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\delta B_p/B_0$",LeftRight,RightTicks);
