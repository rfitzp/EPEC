import graph;
     
size (500, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage5/results.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] r   = A[1];
real[] wnl = A[3];
real[] tt  = A[4]*1.e3;

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
	q4.push (tt[j]);
	i4.push (wnl[j]);
      }
    if ((int) m[j] == 5)
      {
	q5.push (tt[j]);
	i5.push (wnl[j]);
      }
    if ((int) m[j] == 6)
      {
	q6.push (tt[j]);
	i6.push (wnl[j]);
      }
    if ((int) m[j] == 7)
      {
	q7.push (tt[j]);
	i7.push (wnl[j]);
      }
    if ((int) m[j] == 8)
      {
	q8.push (tt[j]);
	i8.push (wnl[j]);
      }
    if ((int) m[j] == 9)
      {
	q9.push (tt[j]);
	i9.push (wnl[j]);
      }
  if ((int) m[j] == 10)
      {
	q10.push (tt[j]);
	i10.push (wnl[j]);
      }
  
    if ((int) m[j] == 11)
      {
	q11.push (tt[j]);
	i11.push (wnl[j]);
      }
    if ((int) m[j] == 12)
      {
	q12.push (tt[j]);
	i12.push (wnl[j]);
      }
    if ((int) m[j] == 13)
      {
	q13.push (tt[j]);
	i13.push (wnl[j]);
      }
    if ((int) m[j] == 14)
      {
	q14.push (tt[j]);
	i14.push (wnl[j]);
      }
    if ((int) m[j] == 15)
      {
	q15.push (tt[j]);
	i15.push (wnl[j]);
      }
    if ((int) m[j] == 16)
      {
	q16.push (tt[j]);
	i16.push (wnl[j]);
      }
  }

pen s;

if (q16.length > 0)
   {
     s = palegreen + solid + 0.5;
     draw (graph (q16, i16), s);
   }
if (q15.length > 0)
   {
    s = grey + solid + 0.5;
    draw (graph (q15, i15), s);
   }	
if (q14.length > 0)
   {
    s = orange + solid + 0.5;
   draw (graph (q14, i14), s);
   }	
if (q13.length > 0)
   {
    s = purple + solid + 0.5;
    draw (graph (q13, i13), s);
   }
if (q12.length > 0)
   {
    s = pink + solid + 0.5;
    draw (graph (q12, i12), s);
   }	
if (q11.length > 0)
   {
    s = brown + solid + 0.5;
    draw (graph (q11, i11), s);
   }	
if (q10.length > 0)
   {
    s = magenta + solid + 0.5;
    draw (graph (q10, i10), s);
   }	
if (q9.length > 0)
   {
    s = cyan + solid + 0.5;
    draw (graph (q9, i9),   s);
   }
if (q8.length > 0)
   {
    s = yellow + solid + 0.5;
    draw (graph (q8, i8),   s);
   }
if (q7.length > 0)
   {
    s = blue + solid + 0.5;
    draw (graph (q7, i7),   s);
   }
if (q6.length > 0)
   {
    s = green + solid + 0.5;
    draw (graph (q6, i6),   s);
   }
if (q5.length > 0)
   {
    s = red + solid + 0.5;
    draw (graph (q5, i5),   s);
   }
if (q4.length > 0)
   {
    s = black + solid + 0.5;
    draw (graph (q4, i4),   s);
   }

s = dotted + black + 1;
ylimits (-50., 50., Crop);
yequals (0., s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t ({\rm ms})$", BottomTop, LeftTicks);
yaxis ("$\varpi({\rm krad/s})$", LeftRight, RightTicks);
