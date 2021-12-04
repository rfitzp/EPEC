import graph;

file    in = input ("velocity.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

real[] x  = A[0];
real[] ve = A[1];
real[] v  = 0.5 * ve;
real[] vi = 0.0 * ve;

size (750, 500, IgnoreAspect);

pen s0;

limits ((-8.,-0.1),(8.,1.7),Crop);

s0 = black + dashed  + 1.;
xequals (2., s0);
xequals (-2., s0);

s0 = red + 1.5;
draw (graph (x, ve), s0);

s0 = black + 1.5;
draw (graph (x, v), s0);

s0 = blue + 1.5;
draw (graph (x, vi), s0);

pen qq = fontsize (20.);
defaultpen (qq);
xaxis ("$X$", BottomTop, LeftTicks);
yaxis ("$\hat{V}_\theta$", LeftRight, RightTicks);


