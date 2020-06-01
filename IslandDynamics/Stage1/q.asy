import graph;
     
size (750, 500, IgnoreAspect);

file    in = input ("q.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] q0 = A[0];
real[] q  = A[1];
real[] qa = A[2];
real[] t  = A[3];

pen s = white + dashed + 1.;

draw (graph (t, q0), s, marker (scale (1.mm)*polygon (3),  red));
draw (graph (t, q),  s, marker (scale (1.mm)*polygon (10), green));
draw (graph (t, qa), s, marker (scale (1.mm)*polygon (4),  blue));

s = dotted + black + 1;
yequals (0.5, s);
yequals (1.0, s);
yequals (1.5, s);
yequals (2.0, s);
yequals (2.5, s);
yequals (3.0, s);
yequals (3.5, s);
yequals (4.0, s);
yequals (4.5, s);
yequals (5.0, s);
yequals (5.5, s);
yequals (6.0, s);
yequals (6.5, s);
yequals (7.0, s);
yequals (7.5, s);
yequals (8.0, s);
yequals (8.5, s);

//ylimits (0.,5., Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",           BottomTop, LeftTicks);
yaxis ("$q_0, q_{95}, q_a$", LeftRight, RightTicks ("$% #.1f$", Step = 1.0, step = 0.1));
