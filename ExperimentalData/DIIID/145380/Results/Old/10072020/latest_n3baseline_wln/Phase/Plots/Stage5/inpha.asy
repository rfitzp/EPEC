import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/mirnov.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p     = A[0];
real[] chi1  = A[1];
real[] chi2  = A[2];
real[] chi3  = A[3];
real[] chi4  = A[4];

real[] pha;
for (int i = 0; i < p.length; i += 1)
  {
  pha.push (atan2 (chi2[i], chi1[i]) /pi);
  }
 
pen s = black + 1.5;	
draw(graph(p,pha),s);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\varphi_{\rm in}/\pi$",LeftRight,RightTicks);
