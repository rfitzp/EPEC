import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage4/ufile.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] q = A[0];
real[] p = A[1];
real[] w = A[2];

file    in1 = input("../../Outputs/Stage4/mfile.txt").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] q1 = A1[0];
real[] p1 = A1[1];
real[] w1 = A1[2];

file    in2 = input("../../Outputs/Stage4/lfile.txt").line();
real[][] A2 = in2.dimension (0,0);
A2          = transpose(A2);
     
real[] q2 = A2[0];
real[] p2 = A2[1];
real[] w2 = A2[2];

pen s = red + solid + 1.5;	
draw(graph(p,w),s,marker(scale(2.0mm)*polygon(4)));

s = black + solid + 1.5;	
draw(graph(p1,w1),s,marker(scale(2.0mm)*polygon(4)));

s = blue + solid + 1.5;	
draw(graph(p2,w2),s,marker(scale(2.0mm)*polygon(4)));

xlimits (0.85,1.0,Crop);
ylimits (0, 0.035,Crop);

s = dotted + black + 2;

xequals (0.945, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$W_{\rm vac}$",LeftRight,RightTicks);
