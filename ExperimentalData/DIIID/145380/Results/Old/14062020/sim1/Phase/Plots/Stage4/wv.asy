import graph;
     
size(500,500,IgnoreAspect);

file    in = input("../../Outputs/Stage4/vac.txt").line();
real[][] A = in.dimension (0,0);
int nres   = A[0].length - 1;
A          = transpose(A);
     
real[] p = A[0];
pen s;

for (int i = 0; i < nres; i += 1)
  {
    if (i == 0)
      s = black;
    else if (i == 1)
      s = red;
    else if (i == 2)
      s = green;
    else if (i == 3)
      s = blue;
    else if (i == 4)
      s = yellow;
    else if (i == 5)
      s = cyan;
    else if (i == 6)
      s = magenta;
    else if (i == 7)
      s = pink;
    else if (i == 8)
      s = brown;
    else if (i == 9)
      s = purple;
    else if (i == 10)
      s = orange;
    draw (graph (p, A[i+1]), s);
  }
xlimits (0., 4.0, Crop);

s = dotted + black + 1;
yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Delta}\phi_{UL}/\pi$",BottomTop,LeftTicks);
yaxis("$W_k({\rm m})$",LeftRight,RightTicks);
