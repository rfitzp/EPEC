import graph;
     
size (750, 500, IgnoreAspect);

file    in = input ("q.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] q0 = A[0];
real[] q  = A[1];
real[] qa = A[2];
real[] ql = A[3];
real[] t  = A[4];

pen s = white + dashed + 1.;

draw (graph (t, q0), s, marker (scale (1.mm)*polygon (3),  red));
draw (graph (t, q),  s, marker (scale (1.mm)*polygon (10), green));
draw (graph (t, ql), s, marker (scale (1.mm)*polygon (4),  blue));
draw (graph (t, qa), s, marker (scale (1.mm)*polygon (4),  black));

s = dotted + black + 1;

yequals (0.00000, s);
yequals (1.00000, s);
yequals (1.33333, s);
yequals (1.66666, s);
yequals (2.00000, s);
yequals (2.33333, s);
yequals (2.66666, s);
yequals (3.00000, s);
yequals (3.33333, s);
yequals (3.66666, s);
yequals (4.00000, s);
yequals (4.33333, s);
yequals (4.66666, s);
yequals (5.00000, s);
yequals (5.33333, s);
yequals (5.66666, s);
yequals (6.00000, s);
yequals (6.33333, s);
yequals (6.66666, s);
yequals (7.00000, s);
yequals (7.33333, s);
yequals (7.66666, s);
yequals (8.00000, s);
yequals (8.33333, s);
yequals (8.66666, s);

//ylimits (0.,5., Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",           BottomTop, LeftTicks);
yaxis ("$q_0, q_{95}, q_l, q_a$", LeftRight, RightTicks ("$% #.1f$", Step = 1.0, step = 0.1));
