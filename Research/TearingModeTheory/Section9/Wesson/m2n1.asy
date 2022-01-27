import graph;

file    in = input ("m2n1r12.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

real[] qa    = A[0];
real[] rs    = A[1];
real[] delta = A[2];
real[] w     = A[3];

file    in1 = input ("m2n1r11.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

real[] qa1    = A1[0];
real[] rs1    = A1[1];
real[] delta1 = A1[2];
real[] w1     = A1[3];

file    in2 = input ("m2n1r10.out").line();
real[][] A2 = in2.dimension (0,0);
A2          = transpose(A2);

real[] qa2    = A2[0];
real[] rs2    = A2[1];
real[] delta2 = A2[2];
real[] w2     = A2[3];

size (750, 500, IgnoreAspect);

pen s0 = red + 1.5;
draw (graph (rs, w), s0);

pen s0 = green + 1.5;
draw (graph (rs1, w1), s0);

pen s0 = blue + 1.5;
draw (graph (rs2, w2), s0);

limits ((0.3,0.),(1.,0.5),Crop);

s0 = black + dotted + 0.5;
yequals (0., s0);

pen qq = fontsize (20.);
defaultpen (qq);
xaxis ("$r_s/a$", BottomTop, LeftTicks);
yaxis ("$W_{pw}/a$", LeftRight, RightTicks);


