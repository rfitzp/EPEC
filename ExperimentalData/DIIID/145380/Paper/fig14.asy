import graph;
     
size (500, 250, IgnoreAspect);

real alpha = 5.20603250270005;
real beta  = -0.000445721197123001;

real f(real x) {return -alpha-beta*x;}

real q0 = f(2450.);
real q1 = f(4351.);

fill((f(2840),0.005)--(f(2980),0.005)--(f(2980),0.995)--(f(2840),0.995)--cycle, paleyellow);
fill((f(3320),0.005)--(f(3560),0.005)--(f(3560),0.995)--(f(3320),0.995)--cycle, paleyellow);
fill((f(3880),0.005)--(f(4200),0.005)--(f(4200),0.995)--(f(3880),0.995)--cycle, paleyellow);

pen p = brown + 6.;
draw ((-3.922,0.75)--(-3.846,0.75),p);
p = magenta + 6.;
draw ((-3.768,0.75)--(-3.680,0.75),p);

p = brown + 6.;
draw ((-3.931,0.5)--(-3.844,0.5),p);
p = magenta + 6.;
draw ((-3.772,0.5)--(-3.643,0.5),p);
p = cyan + 6.;
draw ((-3.455,0.5)--(-3.350,0.5),p);

p = magenta + 6.;
draw ((-3.771,0.25)--(-3.695,0.25),p);
p = cyan + 6.;
draw ((-3.451,0.25)--(-3.381,0.25),p);

pen pp = fontsize(12.);
label(Label("Case 1",(-4.05,0.75)));

pen pp = fontsize(12.);
label(Label("Case 2",(-4.05,0.5)));

pen pp = fontsize(12.);
label(Label("Case 3",(-4.05,0.25)));

limits ((q0,0.),(q1,1.),Crop);

scale(Linear(-1),Linear);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis("$\overline{q_{95}}$", BottomTop, LeftTicks("$% #.1f$",Step=0.1));
yaxis(LeftRight);
