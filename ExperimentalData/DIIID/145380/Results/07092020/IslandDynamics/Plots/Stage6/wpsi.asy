import graph;
     
size (750, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/omega.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] tt  = A[4];
real[] w   = A[7];
real[] psi = A[6];

int N = m.length;
real[] t4,  wl4,  wu4,  wx4,   wy4;
real[] t5,  wl5,  wu5,  wx5,   wy5;
real[] t6,  wl6,  wu6,  wx6,   wy6;
real[] t7,  wl7,  wu7,  wx7,   wy7;
real[] t8,  wl8,  wu8,  wx8,   wy8;
real[] t9,  wl9,  wu9,  wx9,   wy9;
real[] t10, wl10, wu10, wx10,  wy10;
real[] t11, wl11, wu11, wx11,  wy11;
real[] t12, wl12, wu12, wx12,  wy12;
real[] t13, wl13, wu13, wx13,  wy13;
real[] t14, wl14, wu14, wx14,  wy14;
real[] t15, wl15, wu15, wx15,  wy15;
real[] t16, wl16, wu16, wx16,  wy16;
for (int j = 0; j < N; ++j)
  {
   if ((int) m[j] == 4)
      {
	t4.push (tt[j]);
	wl4.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu4.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx4.push (psi[j]-w[j]/2.);
	wy4.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 5)
      {
	t5.push (tt[j]);
	wl5.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu5.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
        wx5.push (psi[j]-w[j]/2.);
	wy5.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 6)
      {
	t6.push (tt[j]);
	wl6.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu6.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
        wx6.push (psi[j]-w[j]/2.);
	wy6.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 7)
      {
	t7.push (tt[j]);
	wl7.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu7.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
        wx7.push (psi[j]-w[j]/2.);
	wy7.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 8)
      {
	t8.push (tt[j]);
	wl8.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu8.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
        wx8.push (psi[j]-w[j]/2.);
	wy8.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 9)
      {
	t9.push (tt[j]);
	wl9.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu9.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
        wx9.push (psi[j]-w[j]/2.);
	wy9.push (psi[j]+w[j]/2.);	
      }
  if ((int) m[j] == 10)
      {
	t10.push (tt[j]);
	wl10.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu10.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
        wx10.push (psi[j]-w[j]/2.);
	wy10.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 11)
      {
	t11.push (tt[j]);
 	wl11.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu11.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
        wx11.push (psi[j]-w[j]/2.);
	wy11.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 12)
      {
	t12.push (tt[j]);
  	wl12.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu12.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx12.push (psi[j]-w[j]/2.);
	wy12.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 13)
      {
	t13.push (tt[j]);
 	wl13.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu13.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx13.push (psi[j]-w[j]/2.);
	wy13.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 14)
      {
	t14.push (tt[j]);
     	wl14.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu14.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx14.push (psi[j]-w[j]/2.);
	wy14.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 15)
      {
	t15.push (tt[j]);
	wl15.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu15.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx15.push (psi[j]-w[j]/2.);
	wy15.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 16)
      {
	t16.push (tt[j]);
	wl16.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu16.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx16.push (psi[j]-w[j]/2.);
	wy16.push (psi[j]+w[j]/2.);	
      }
  }

pen s;

fill((2840,0.8505)--(2980,0.8505)--(2980,0.9995)--(2840,0.9995)--cycle, paleyellow);
fill((3320,0.8505)--(3560,0.8505)--(3560,0.9995)--(3320,0.9995)--cycle, paleyellow);
fill((3880,0.8505)--(4200,0.8505)--(4200,0.9995)--(3880,0.9995)--cycle, paleyellow);

if (t16.length > 0)
   {
     s = palegreen + solid + 0.5;
     //draw (graph (t16, wl16), s);
     //draw (graph (t16, wu16), s);
     //draw (graph (t16, wx16), s);
     //draw (graph (t16, wy16), s);
   }
if (t15.length > 0)
   {
     s = grey + solid + 0.5;
    // draw (graph (t15, wl15), s);
     //draw (graph (t15, wu15), s);
     //draw (graph (t15, wx15), s);
    // draw (graph (t15, wy15), s);
    }	
if (t11.length > 0)
   {
     s = brown + solid + 0.5;
     draw (graph (t11, wl11), s);
     draw (graph (t11, wu11), s);
     draw (graph (t11, wx11), s);
     draw (graph (t11, wy11), s);
   }	
if (t12.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t12, wl12), s);
     draw (graph (t12, wu12), s);
     draw (graph (t12, wx12), s);
     draw (graph (t12, wy12), s);
  } 	
if (t13.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t13, wl13), s);
     draw (graph (t13, wu13), s);
     draw (graph (t13, wx13), s);
     draw (graph (t13, wy13), s);
   }
if (t14.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t14, wl14), s);
     draw (graph (t14, wu14), s);
     draw (graph (t14, wx14), s);
     draw (graph (t14, wy14), s);
   }	
if (t10.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t10, wl10), s);
     draw (graph (t10, wu10), s);
     draw (graph (t10, wx10), s);
     draw (graph (t10, wy10), s);
   }	
if (t9.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t9, wl9), s);
     draw (graph (t9, wu9), s);
     draw (graph (t9, wx9), s);
     draw (graph (t9, wy9), s);
   }
if (t8.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t8, wl8), s);
     draw (graph (t8, wu8), s);
     draw (graph (t8, wx8), s);
     draw (graph (t8, wy8), s);
   }
if (t7.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t7, wl7), s);
     draw (graph (t7, wu7), s);
     draw (graph (t7, wx7), s);
     draw (graph (t7, wy7), s);
   }
if (t6.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t6, wl6), s);
     draw (graph (t6, wu6), s);
     draw (graph (t6, wx6), s);
     draw (graph (t6, wy6), s);
   }
if (t5.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t5, wl5), s);
     draw (graph (t5, wu5), s);
     draw (graph (t5, wx5), s);
     draw (graph (t5, wy5), s);
   }
if (t4.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t4, wl4), s);
     draw (graph (t4, wu4), s);
     draw (graph (t4, wx4), s);
     draw (graph (t4, wy4), s);
   }

s = dotted + black + 1;
ylimits (0.85, 1., Crop);

s = dotted + black + 2;
yequals (0.925, s);

//yequals (0., s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t ({\rm ms})$", BottomTop, LeftTicks);
yaxis ("${\mit\Psi}_N$", LeftRight, RightTicks);
