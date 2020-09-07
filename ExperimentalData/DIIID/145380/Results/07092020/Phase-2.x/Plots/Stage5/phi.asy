import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/phi.txt").line();
real[][] A = in.dimension (0,0);
int nres   = A[0].length - 1;
A          = transpose(A);
     
real[] p = A[0];
pen s;

for (int i = 0; i < nres; i += 1)
  {
    if (i%4 == 0)
      s = black;
      else if (i%4 == 1)
	s = red;
      else if (i%4 == 2)
	s = green;
      else if (i%4 == 3)
	s = blue;
    draw (graph (p, A[i+1]), s);
  }
s = dotted + black + 1;
yequals (0., s);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\varphi_k - \zeta_k$",LeftRight,RightTicks);
