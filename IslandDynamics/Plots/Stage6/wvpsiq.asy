import graph;
     
size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/omega.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] tt  = A[4];
real[] w   = A[8];
real[] psi = A[6];
real[] q   = A[11];

int N = m.length;
real[] q2, t2,  wl2,  wu2;
real[] q3, t3,  wl3,  wu3;
real[] q4, t4,  wl4,  wu4;
real[] q5, t5,  wl5,  wu5;
real[] q6, t6,  wl6,  wu6;
real[] q7, t7,  wl7,  wu7;
real[] q8, t8,  wl8,  wu8;
real[] q9, t9,  wl9,  wu9;
real[] q10, t10, wl10, wu10;
real[] q11, t11, wl11, wu11;
real[] q12, t12, wl12, wu12;
real[] q13, t13, wl13, wu13;
real[] q14, t14, wl14, wu14;
real[] q15, t15, wl15, wu15;
real[] q16, t16, wl16, wu16;
for (int j = 0; j < N; ++j)
  {
   if ((int) m[j] == 2)
      {
        q2.push (q[j]);
	t2.push (tt[j]);
	wl2.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu2.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 3)
      {
        q3.push (q[j]);
	t3.push (tt[j]);
	wl3.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu3.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
        if ((int) m[j] == 4)
      {
        q4.push (q[j]);
	t4.push (tt[j]);
	wl4.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu4.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 5)
      {
        q5.push (q[j]);
	t5.push (tt[j]);
	wl5.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu5.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 6)
      {
        q6.push (q[j]);
	t6.push (tt[j]);
	wl6.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu6.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 7)
      {
        q7.push (q[j]);
	t7.push (tt[j]);
	wl7.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu7.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 8)
      {
        q8.push (q[j]);
	t8.push (tt[j]);
	wl8.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu8.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 9)
      {
        q9.push (q[j]);
	t9.push (tt[j]);
	wl9.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu9.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
  if ((int) m[j] == 10)
      {
        q10.push (q[j]);
	t10.push (tt[j]);
	wl10.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu10.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 11)
      {
        q11.push (q[j]);
	t11.push (tt[j]);
 	wl11.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu11.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
      }
    if ((int) m[j] == 12)
      {
        q12.push (q[j]);
	t12.push (tt[j]);
  	wl12.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu12.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        }
    if ((int) m[j] == 13)
      {
        q13.push (q[j]);
	t13.push (tt[j]);
 	wl13.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu13.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
	}
    if ((int) m[j] == 14)
      {
        q14.push (q[j]);
	t14.push (tt[j]);
     	wl14.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu14.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
	}
    if ((int) m[j] == 15)
      {
        q15.push (q[j]);
	t15.push (tt[j]);
	wl15.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu15.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
      }
    if ((int) m[j] == 16)
      {
        q16.push (q[j]);
	t16.push (tt[j]);
	wl16.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu16.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
      }
  }

pen s;

if (t16.length > 0)
   {
     s = palegreen + solid + 0.5;
     draw (graph (q16, wl16), s);
     draw (graph (q16, wu16), s);
   }
if (t15.length > 0)
   {
     s = grey + solid + 0.5;
     draw (graph (q15, wl15), s);
     draw (graph (q15, wu15), s);
    }	
if (t14.length > 0)
   {
     s = orange + solid + 0.5;
    draw (graph (q14, wl14), s);
    draw (graph (q14, wu14), s);
   }	
if (t13.length > 0)
   {  
     s = purple + solid + 0.5;
    draw (graph (q13, wl13), s);
    draw (graph (q13, wu13), s);
   }
if (t12.length > 0)
   {
     s = pink + solid + 0.5;
    draw (graph (q12, wl12), s);
    draw (graph (q12, wu12), s);
  } 	
if (t11.length > 0)
   {
     s = brown + solid + 0.5;
    draw (graph (q11, wl11), s);
    draw (graph (q11, wu11), s);
   }	
if (t10.length > 0)
   {
     s = magenta + solid + 0.5;
    draw (graph (q10, wl10), s);
    draw (graph (q10, wu10), s);
   }	
if (t9.length > 0)
   {
     s = cyan + solid + 0.5;
    draw (graph (q9, wl9), s);
    draw (graph (q9, wu9), s);
   }
if (t8.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (q8, wl8), s);
     draw (graph (q8, wu8), s);
   }
if (t7.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (q7, wl7), s);
     draw (graph (q7, wu7), s);
   }
if (t6.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (q6, wl6), s);
     draw (graph (q6, wu6), s);
   }
if (t5.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (q5, wl5), s);
     draw (graph (q5, wu5), s);
   }
if (t4.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (q4, wl4), s);
     draw (graph (q4, wu4), s);
   }
if (t3.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (q3, wl3), s);
     draw (graph (q3, wu3), s);
   }
if (t2.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (q2, wl2), s);
     draw (graph (q2, wu2), s);
   }

s = dotted + black + 1;
ylimits (0.8, 1., Crop);
//yequals (0., s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$q_{95}$", BottomTop, LeftTicks);
yaxis ("${\mit\Psi}_N$", LeftRight, RightTicks);
