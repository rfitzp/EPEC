import graph;
     
size (750, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/omega.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] tt  = A[4];
real[] w   = A[8];
real[] psi = A[6];

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
	t4.push (tt[j]);
	wl4.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu4.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 5)
      {
	t5.push (tt[j]);
	wl5.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu5.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 6)
      {
	t6.push (tt[j]);
	wl6.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu6.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 7)
      {
	t7.push (tt[j]);
	wl7.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu7.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 8)
      {
	t8.push (tt[j]);
	wl8.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu8.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 9)
      {
	t9.push (tt[j]);
	wl9.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu9.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
  if ((int) m[j] == 10)
      {
	t10.push (tt[j]);
	wl10.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu10.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);	
      }
    if ((int) m[j] == 11)
      {
	t11.push (tt[j]);
 	wl11.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu11.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
      }
    if ((int) m[j] == 12)
      {
	t12.push (tt[j]);
  	wl12.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu12.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        }
    if ((int) m[j] == 13)
      {
	t13.push (tt[j]);
 	wl13.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu13.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
	}
    if ((int) m[j] == 14)
      {
	t14.push (tt[j]);
     	wl14.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu14.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
	}
    if ((int) m[j] == 15)
      {
	t15.push (tt[j]);
	wl15.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu15.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
      }
    if ((int) m[j] == 16)
      {
	t16.push (tt[j]);
	wl16.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu16.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
      }
  }

pen s;

if (t16.length > 0)
   {
     s = palegreen + solid + 0.5;
     draw (graph (t16, wl16), s);
     draw (graph (t16, wu16), s);
   }
if (t15.length > 0)
   {
     s = grey + solid + 0.5;
     draw (graph (t15, wl15), s);
     draw (graph (t15, wu15), s);
    }	
if (t14.length > 0)
   {
     s = orange + solid + 0.5;
    draw (graph (t14, wl14), s);
    draw (graph (t14, wu14), s);
   }	
if (t13.length > 0)
   {  
     s = purple + solid + 0.5;
    draw (graph (t13, wl13), s);
    draw (graph (t13, wu13), s);
   }
if (t12.length > 0)
   {
     s = pink + solid + 0.5;
    draw (graph (t12, wl12), s);
    draw (graph (t12, wu12), s);
  } 	
if (t11.length > 0)
   {
     s = brown + solid + 0.5;
    draw (graph (t11, wl11), s);
    draw (graph (t11, wu11), s);
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
     draw (graph (t5, wl5), s);
     draw (graph (t5, wu5), s);
   }
if (t4.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t4, wl4), s);
     draw (graph (t4, wu4), s);
   }

s = dotted + black + 1;
ylimits (0.8, 1., Crop);
//yequals (0., s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t ({\rm ms})$", BottomTop, LeftTicks);
yaxis ("${\mit\Psi}_N$", LeftRight, RightTicks);