import graph;
     
size (750, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/q.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] q0 = A[0];
real[] q  = A[1];
real[] qa = A[2];
real[] ql = A[3];
real[] t  = A[4];

fill((2840,0.05)--(2980,0.05)--(2980,8.95)--(2840,8.95)--cycle, paleyellow);
fill((3320,0.05)--(3560,0.05)--(3560,8.95)--(3320,8.95)--cycle, paleyellow);
fill((3880,0.05)--(4200,0.05)--(4200,8.95)--(3880,8.95)--cycle, paleyellow);

pen s = white + dashed + 1.;

//draw (graph (t, q0), s, marker (scale (0.5mm)*polygon (3),  red));
draw (graph (t, q),  s, marker (scale (0.5mm)*polygon (10), green));
//draw (graph (t, ql), s, marker (scale (0.5mm)*polygon (4),  blue));
//draw (graph (t, qa), s, marker (scale (0.5mm)*polygon (4),  black));

s = dotted + black + 1;

//yequals (0.00000, s);
//yequals (0.33333, s);
//yequals (0.66666, s);
//yequals (1.00000, s);
//yequals (1.33333, s);
//yequals (1.66666, s);
//yequals (2.00000, s);
//yequals (2.33333, s);
//yequals (2.66666, s);
//yequals (3.00000, s);
//yequals (3.33333, s);
//yequals (3.66666, s);
//yequals (4.00000, s);
//yequals (4.33333, s);
//yequals (4.66666, s);
//yequals (5.00000, s);
//yequals (5.33333, s);
//yequals (5.66666, s);
//yequals (6.00000, s);
//yequals (6.33333, s);
//yequals (6.66666, s);
//yequals (7.00000, s);
//yequals (7.33333, s);
//yequals (7.66666, s);
//yequals (8.00000, s);
//yequals (8.33333, s);
//yequals (8.66666, s);

yequals (3.95,s);
yequals (3.6,s);
yequals (3.48,s);

limits ((2400.,3.2),(4300.,4.2), Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t({\rm ms})$",           BottomTop, LeftTicks);
yaxis ("$q$", LeftRight, RightTicks ("$% #.1f$", Step = 1.0, step = 0.1));
