import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage5/omega.txt").line();
real[][] A = in.dimension (0,0);
int nres   = A[0].length - 1;
A          = transpose(A);
     
real[] p = A[0];
pen s;

for (int i = 0; i < nres; i += 1)
  {
    if (i%11 == 0)
      s = black;
    else if (i%11 == 1)
      s = red;
    else if (i%11 == 2)
      s = green;
    else if (i%11 == 3)
      s = blue;
    else if (i%11 == 4)
      s = yellow;
    else if (i%11 == 5)
      s = cyan;
    else if (i%11 == 6)
      s = magenta;
    else if (i%11 == 7)
      s = brown;
    else if (i%11 == 8)
      s = pink;
    else if (i%11 == 9)
      s = purple;
    else if (i%11 == 10)
      s = orange;
      
    draw (graph (p, A[i+1]), s);
  }

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$t ({\rm s})$",BottomTop,LeftTicks);
yaxis("$\omega_k ({\rm krad/s})$",LeftRight,RightTicks);
