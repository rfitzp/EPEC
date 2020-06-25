import graph;
     
size(750,500,IgnoreAspect);

file    in = input("../../Outputs/Stage3/profiles.txt").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] psi = A[0];
real[] r   = A[1];
real[] ne  = A[2];
real[] te  = A[4];
real[] ti  = A[8];
real[] nii = A[10]*10.;

pen s = green + solid + 1.5;	
draw(graph(psi,ne),s);
s = red + solid + 1.5;	
draw(graph(psi,te),s);
s = blue + solid + 1.5;	
draw(graph(psi,ti),s);
s = cyan + solid + 1.5;	
draw(graph(psi,nii),s);

s = dotted + black + 1;

limits ((0.,0.),(1.,6.5),Crop);

//yequals (0., s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("${\mit\Psi}_N$",BottomTop,LeftTicks);
qq = fontsize(15.);
yaxis("$\textcolor{green}{n_e}(10^{19}\,{\rm m}^{-3}), \textcolor{cyan}{n_I}(10^{18}\,{\rm m}^{-3}), \textcolor{red}{T_e}({\rm keV}), \textcolor{blue}{T_i}({\rm keV})$",LeftRight,RightTicks);
