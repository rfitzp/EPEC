import graph;
     
size(1000,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] ne  = A[2];
real[] te  = A[4];
real[] ni  = A[6];
real[] ti  = A[8];
real[] nii = A[10]*10.;
real[] nn  = A[21]*100.;

pen s = red + solid + 4.;	
draw(graph(psi,ne),s);

draw ((0.7,6.)--(0.8,6.),s);

s = green + solid + 4.0;	
draw(graph(psi,te),s);

draw ((0.7,5.5)--(0.8,5.5),s);

s = blue + solid + 4.0;	
draw(graph(psi,ti),s);

draw ((0.7,5.)--(0.8,5.),s);

s = cyan + solid + 4.0;	
draw(graph(psi,nii),s);

draw ((0.7,4.5)--(0.8,4.5),s);

s = magenta + solid + 3.0;	
draw(graph(psi,nn),s);

draw ((0.7,4.)--(0.8,4.),s);

s = dotted + black + 1;

limits ((0.,0.),(1.,6.5),Crop);

//yequals (0., s);

s = dotted + black + 2;

xequals (0.945, s);

pen pp = fontsize(30.);
label(Label("$n_e$",(0.83,6.)));
label(Label("$T_e$",(0.83,5.5)));
label(Label("$T_i$",(0.83,5.)));
label(Label("$10\,n_I$",(0.845,4.5)));
label(Label("$100\,\langle n_n\rangle$",(0.86,4.)));

pen qq = fontsize(50.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
yaxis("$n_{e,I,n}$, $T_{e,i}$",LeftRight,RightTicks);
