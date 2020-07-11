import graph;
     
size (500, 500, IgnoreAspect);

file    in = input ("scan.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);

real ILIM = 4.;

ILIM = getreal ("Ilim (kA) ? ");

real[] q = A[0];
real[] m = A[1];
real[] r = A[2];
real[] t = A[3];
real[] i = A[4];
real[] p = A[5];
real[] tt = A[12];

int N = q.length;
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
   if ((int) m[j] == 4 && i[j] < ILIM)
      {
	q4.push (tt[j]);
	i4.push (r[j]);
      }
    if ((int) m[j] == 5 && i[j] < ILIM)
      {
	q5.push (tt[j]);
	i5.push (r[j]);
      }
    if ((int) m[j] == 6 && i[j] < ILIM)
      {
	q6.push (tt[j]);
	i6.push (r[j]);
      }
    if ((int) m[j] == 7 && i[j] < ILIM)
      {
	q7.push (tt[j]);
	i7.push (r[j]);
      }
    if ((int) m[j] == 8 && i[j] < ILIM)
      {
	q8.push (tt[j]);
	i8.push (r[j]);
      }
    if ((int) m[j] == 9 && i[j] < ILIM)
      {
	q9.push (tt[j]);
	i9.push (r[j]);
      }
    if ((int) m[j] == 10 && i[j] < ILIM)
      {
	q10.push (tt[j]);
	i10.push (r[j]);
      }
    if ((int) m[j] == 11 && i[j] < ILIM)
      {
	q11.push (tt[j]);
	i11.push (r[j]);
      }
    if ((int) m[j] == 12 && i[j] < ILIM)
      {
	q12.push (tt[j]);
	i12.push (r[j]);
      }
    if ((int) m[j] == 13 && i[j] < ILIM)
      {
	q13.push (tt[j]);
	i13.push (r[j]);
      }
    if ((int) m[j] == 14 && i[j] < ILIM)
      {
	q14.push (tt[j]);
	i14.push (r[j]);
      }
    if ((int) m[j] == 15 && i[j] < ILIM)
      {
	q15.push (tt[j]);
	i15.push (r[j]);
      }
    if ((int) m[j] == 16 && i[j] < ILIM)
      {
	q16.push (tt[j]);
	i16.push (r[j]);
      }
  }

pen s;

s = white + dotted + 0.;
if (q4.length > 0)
   {  
     draw (graph (q4, i4),   s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, black));
   }
if (q5.length > 0)
   {  
     draw (graph (q5, i5),   s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, red));
   }
if (q6.length > 0)
   {  
     draw (graph (q6, i6),   s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, green));
   }
if (q7.length > 0)
   {  
     draw (graph (q7, i7),   s, marker (scale(2.0mm)*polygon(3),   filltype=Fill, blue));
   }
if (q8.length > 0)
   {  
     draw (graph (q8, i8),   s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, yellow));
   }
if (q9.length > 0)
   {  
     draw (graph (q9, i9),   s, marker (scale(2.0mm)*polygon(3), filltype=Fill, cyan));
   }
if (q10.length > 0)
   {  
     draw (graph (q10, i10), s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, magenta));
   }	
if (q11.length > 0)
   {  
     draw (graph (q11, i11), s, marker (scale(2.0mm)*polygon(3), filltype=Fill, brown));
   }	
if (q12.length > 0)
   {  
     draw (graph (q12, i12), s, marker (scale(2.0mm)*polygon(3), filltype=Fill, pink));
   }	
if (q13.length > 0)
   {  
     draw (graph (q13, i13), s, marker (scale(2.0mm)*polygon(3), filltype=Fill, purple));
   }
if (q14.length > 0)
   {  
     draw (graph (q14, i14), s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, orange));
   }	
if (q15.length > 0)
   {  
     draw (graph (q15, i15), s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, grey));
   }	
if (q16.length > 0)
   {  
     draw (graph (q16, i16), s, marker (scale(2.0mm)*polygon(3),  filltype=Fill, grey));
   }	

s = dotted + black + 1;
ylimits (0.7, 1., Crop);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t ({\rm ms})$",                BottomTop, LeftTicks);
yaxis ("$r/a$", LeftRight, RightTicks);