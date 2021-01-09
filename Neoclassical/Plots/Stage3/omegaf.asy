import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/omega.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r    = A[4];
real[] wlin = A[1];
real[] wnl  = A[2];
real[] wEB  = A[3];

pen s = red + solid + 1.5;	
draw(graph(r,wlin),s,marker(scale(2.0mm)*polygon(4)));
s = blue + solid + 1.5;	
draw(graph(r,wnl),s,marker(scale(2.0mm)*polygon(4)));
s = black + solid + 1.5;	
draw(graph(r,wEB),s,marker(scale(2.0mm)*polygon(4)));

xlimits (0.85,1.0,Crop);
ylimits (-150, 550,Crop);

s = dotted + black + 1;
yequals (0., s);

s = dotted + black + 2;

xequals (0.945, s);


pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$\varpi ({\rm krad/s})$",LeftRight,RightTicks);
