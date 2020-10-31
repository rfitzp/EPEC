import graph;
     
size (750, 510, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/omega.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] tt  = A[4];
real[] w   = A[8];
real[] psi = A[6];

real alpha = 5.20603250270005;
real beta  = -0.000445721197123001;

real f(real x) {return -alpha-beta*x;}

int N = m.length;
real[] t4,  wl4,  wu4;
real[] t5,  wl5,  wu5;
real[] t6,  wl6,  wu6;
real[] t7,  wl7,  wu7;
real[] t8,  wl8,  wu8;
real[] t9,  wl9,  wu9;
real[] t10, wl10, wu10;
real[] t11, wl11, wu11;
real[] t12, wl12, wu12;
real[] t13, wl13, wu13;
real[] t14, wl14, wu14;
real[] t15, wl15, wu15;
real[] t16, wl16, wu16;
for (int j = 0; j < N; ++j)
  {
   if ((int) m[j] == 4)
      {
	t4.push (f(tt[j]));
	wl4.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu4.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 5)
      {
	t5.push (f(tt[j]));
	wl5.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu5.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 6)
      {
	t6.push (f(tt[j]));
	wl6.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu6.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 7)
      {
	t7.push (f(tt[j]));
	wl7.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu7.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 8)
      {
	t8.push (f(tt[j]));
	wl8.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu8.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 9)
      {
	t9.push (f(tt[j]));
	wl9.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu9.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
  if ((int) m[j] == 10)
      {
	t10.push (f(tt[j]));
	wl10.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu10.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 11)
      {
	t11.push (f(tt[j]));
 	wl11.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu11.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
      }
    if ((int) m[j] == 12)
      {
	t12.push (f(tt[j]));
  	wl12.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu12.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        }
    if ((int) m[j] == 13)
      {
	t13.push (f(tt[j]));
 	wl13.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu13.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
	}
    if ((int) m[j] == 14)
      {
	t14.push (f(tt[j]));
     	wl14.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu14.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
	}
    if ((int) m[j] == 15)
      {
	t15.push (f(tt[j]));
	wl15.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu15.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
      }
    if ((int) m[j] == 16)
      {
	t16.push (f(tt[j]));
	wl16.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu16.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
      }
  }

pen s;

fill((f(2840),0.8505)--(f(2980),0.8505)--(f(2980),0.9995)--(f(2840),0.9995)--cycle, paleyellow);
fill((f(3320),0.8505)--(f(3560),0.8505)--(f(3560),0.9995)--(f(3320),0.9995)--cycle, paleyellow);
fill((f(3880),0.8505)--(f(4200),0.8505)--(f(4200),0.9995)--(f(3880),0.9995)--cycle, paleyellow);

if (t16.length > 0)
   {
     s = palegreen + solid + 0.5;
     //draw (graph (t16, wl16), s);
    // draw (graph (t16, wu16), s);
   }
if (t15.length > 0)
   {
     s = grey + solid + 0.5;
    // draw (graph (t15, wl15), s);
    // draw (graph (t15, wu15), s);
    }	
if (t11.length > 0)
   {
     s = brown + solid + 0.5;
    draw (graph (t11, wl11), s);
    draw (graph (t11, wu11), s);
   }
   if (t12.length > 0)
   {
     s = pink + solid + 0.5;
    draw (graph (t12, wl12), s);
    draw (graph (t12, wu12), s);
  }
  if (t13.length > 0)
   {
     s = purple + solid + 0.5;
    draw (graph (t13, wl13), s);
    draw (graph (t13, wu13), s);
   }	
if (t14.length > 0)
   {
      s = orange + solid + 0.5;
    draw (graph (t14, wl14), s);
    draw (graph (t14, wu14), s);
     }
if (t10.length > 0)
   {
     s = magenta + solid + 0.5;
    draw (graph (t10, wl10), s);
    draw (graph (t10, wu10), s);
   }	
if (t9.length > 0)
   {
     s = cyan + solid + 0.5;
    draw (graph (t9, wl9), s);
    draw (graph (t9, wu9), s);
   }
if (t8.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t8, wl8), s);
     draw (graph (t8, wu8), s);
   }
if (t7.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t7, wl7), s);
     draw (graph (t7, wu7), s);
   }
if (t6.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t6, wl6), s);
     draw (graph (t6, wu6), s);
   }
if (t5.length > 0)
   {
     s = red + solid + 0.5;
     //draw (graph (t5, wl5), s);
     //draw (graph (t5, wu5), s);
   }
if (t4.length > 0)
   {
     s = black + solid + 0.5;
    // draw (graph (t4, wl4), s);
    // draw (graph (t4, wu4), s);
   }

s = dotted + black + 1;
ylimits (0.85, 1., Crop);
//yequals (0., s);

s = dotted + black + 2;
yequals (0.925, s);

scale(Linear(-1),Linear);

pen qq = fontsize (50.);
defaultpen (qq);
xaxis ("$\overline{q_{95}}$", BottomTop, LeftTicks("$% #.1f$",Step=0.1));
yaxis ("${\mit\Psi}_N$", LeftRight, RightTicks);
      
pen pp = fontsize(40.);
defaultpen (pp);
label(Label("m = "),  (-3.98-0.05+0.16, 1.006),  black);
label(Label("8,"),    (-3.85-0.05+0.1-0.03+0.01, 1.006),  yellow);
label(Label("9,"),    (-3.82-0.05+0.1-0.03+0.02, 1.006),  cyan);
label(Label("10,"),   (-3.79-0.05+0.1-0.03+0.03, 1.006),  magenta);
label(Label("11,"),   (-3.75-0.05+0.1-0.03+0.04, 1.006),  brown);
label(Label("12,"),   (-3.71-0.05+0.1-0.03+0.05, 1.006),  pink);
label(Label("13,"),   (-3.67-0.05+0.1-0.03+0.06, 1.006),  purple);
label(Label("14\phantom{,}"),    (-3.63-0.05+0.1-0.03+0.07, 1.0056),  orange);