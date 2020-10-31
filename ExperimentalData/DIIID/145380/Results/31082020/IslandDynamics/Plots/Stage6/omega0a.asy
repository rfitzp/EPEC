import graph;
     
size (750, 510, IgnoreAspect);

file    in = input ("../../Outputs/Stage6/omega0.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] r   = A[1];
real[] wl  = A[2];
real[] wnl = A[3];
real[] web = A[4];
real[] tt  = A[5];

real alpha = 5.20603250270005;
real beta  = -0.000445721197123001;

real f(real x) {return -alpha-beta*x;}

int N = m.length;
real[] q4,  i4;
real[] q5,  i5;
real[] q6,  i6;
real[] q7,  i7;
real[] q8,  i8;
real[] q9,  i9;
real[] q10, i10;
real[] q11, i11;
real[] q12, i12;
real[] q13, i13;
real[] q14, i14;
real[] q15, i15;
real[] q16, i16;
for (int j = 0; j < N; ++j)
  {
   if ((int) m[j] == 4)
      {
	q4.push (f(tt[j]));
	i4.push (wnl[j]);
      }
    if ((int) m[j] == 5)
      {
	q5.push (f(tt[j]));
	i5.push (wnl[j]);
      }
    if ((int) m[j] == 6)
      {
	q6.push (f(tt[j]));
	i6.push (wnl[j]);
      }
    if ((int) m[j] == 7)
      {
	q7.push (f(tt[j]));
	i7.push (wnl[j]);
      }
    if ((int) m[j] == 8)
      {
	q8.push (f(tt[j]));
	i8.push (wnl[j]);
      }
    if ((int) m[j] == 9)
      {
	q9.push (f(tt[j]));
	i9.push (wnl[j]);
      }
   if ((int) m[j] == 10)
      {
	q10.push (f(tt[j]));
	i10.push (wnl[j]);
      }
    if ((int) m[j] == 11)
      {
	q11.push (f(tt[j]));
	i11.push (wnl[j]);
      }
    if ((int) m[j] == 12)
      {
	q12.push (f(tt[j]));
	i12.push (wnl[j]);
      }
    if ((int) m[j] == 13)
      {
	q13.push (f(tt[j]));
	i13.push (wnl[j]);
      }
    if ((int) m[j] == 14)
      {
	q14.push (f(tt[j]));
	i14.push (wnl[j]);
      }
    if ((int) m[j] == 15)
      {
	q15.push (f(tt[j]));
	i15.push (wnl[j]);
      }
    if ((int) m[j] == 16)
      {
	q16.push (f(tt[j]));
	i16.push (wnl[j]);
      }
  }

pen s;

fill((f(2840),-49.5)--(f(2980),-49.5)--(f(2980),49.5)--(f(2840),49.5)--cycle, paleyellow);
fill((f(3320),-49.5)--(f(3560),-49.5)--(f(3560),49.5)--(f(3320),49.5)--cycle, paleyellow);
fill((f(3880),-49.5)--(f(4200),-49.5)--(f(4200),49.5)--(f(3880),49.5)--cycle, paleyellow);

s = white + dotted + 0.5;
if (q16.length > 0)
   {  
 //    draw (graph (q16, i16), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, palegreen));
   }	
if (q15.length > 0)
   {  
 //   draw (graph (q15, i15), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, grey));
   }	
if (q14.length > 0)
   {  
   draw (graph (q14, i14), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, orange));
   }	
if (q13.length > 0)
   {  
    draw (graph (q13, i13), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, purple));
   }
if (q12.length > 0)
   {  
    draw (graph (q12, i12), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, pink));
  }	
if (q11.length > 0)
   {  
    draw (graph (q11, i11), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, brown));
   }	
if (q10.length > 0)
   {  
   draw (graph (q10, i10), s, marker (scale(0.5mm)*polygon(3), filltype=Fill, magenta));
   }	
if (q9.length > 0)
   {  
     draw (graph (q9, i9),   s, marker (scale(0.5mm)*polygon(3),  filltype=Fill, cyan));
   }
if (q8.length > 0)
   {  
     draw (graph (q8, i8),   s, marker (scale(0.5mm)*polygon(3),  filltype=Fill, yellow));
   }
if (q7.length > 0)
   {  
     draw (graph (q7, i7),   s, marker (scale(0.5mm)*polygon(3),  filltype=Fill, blue));
   }
if (q6.length > 0)
   {  
     draw (graph (q6, i6),   s, marker (scale(0.5mm)*polygon(3),  filltype=Fill, green));
   }
if (q5.length > 0)
   {  
     draw (graph (q5, i5),   s, marker (scale(0.5mm)*polygon(3),  filltype=Fill, red));
   }
if (q4.length > 0)
   {  
     draw (graph (q4, i4),   s, marker (scale(0.5mm)*polygon(3),  filltype=Fill, black));
   }

s = dotted + black + 1;
ylimits (-50., 50., Crop);
yequals (0., s);

scale(Linear(-1),Linear);

pen qq = fontsize (50.);
defaultpen (qq);
xaxis ("$\overline{q_{95}}$", BottomTop, LeftTicks("$% #.1f$",Step=0.1));
yaxis ("$\varpi_0({\rm krad/s})$", LeftRight, RightTicks);

pen pp = fontsize(40.);
defaultpen (pp);
label(Label("m = "),  (-3.98-0.05+0.1, 55),  black);
label(Label("5,"),    (-3.94-0.05+0.1+0.01, 55),  red);
label(Label("6,"),    (-3.91-0.05+0.1+0.02, 55),  green);
label(Label("7,"),    (-3.88-0.05+0.1+0.03, 55),  blue);
label(Label("8,"),    (-3.85-0.05+0.1+0.04, 55),  yellow);
label(Label("9,"),    (-3.82-0.05+0.1+0.05, 55),  cyan);
label(Label("10,"),   (-3.79-0.05+0.1+0.06, 55),  magenta);
label(Label("11,"),   (-3.75-0.05+0.1+0.07, 55),  brown);
label(Label("12,"),   (-3.71-0.05+0.1+0.08, 55),  pink);
label(Label("13,"),   (-3.67-0.05+0.1+0.09, 55),  purple);
label(Label("14\phantom{,}"),    (-3.63-0.05+0.1+0.1, 55),  orange);

