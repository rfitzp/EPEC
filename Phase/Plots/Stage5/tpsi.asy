import graph;
     
size (1000, 500, IgnoreAspect);

file    in = input ("../../Outputs/Stage5/results.txt").line();
real[][] A = in.dimension (0, 0);
A          = transpose (A);
     
real[] m   = A[0];
real[] tt  = A[4]*1.e3;
real[] w   = A[10];
real[] psi = A[6];

file    inx = input ("../../../Flux/Outputs/Stage1/Psilim.txt").line();
real[][] Ax = inx.dimension (0, 0);
real[] ppp  = Ax[0];
real psilim = ppp[0];
real psiped = ppp[1];

int N = m.length;
real[] t1,  wl1,  wu1,  wx1,   wy1;
real[] t2,  wl2,  wu2,  wx2,   wy2;
real[] t3,  wl3,  wu3,  wx3,   wy3;
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
real[] t17, wl17, wu17, wx17,  wy17;
real[] t18, wl18, wu18, wx18,  wy18;
real[] t19, wl19, wu19, wx19,  wy19;

real[] t20, wl20, wu20, wx20,  wy20;
real[] t20, wl20, wu20, wx20,  wy20;
real[] t21, wl21, wu21, wx21,  wy21;
real[] t22, wl22, wu22, wx22,  wy22;
real[] t23, wl23, wu23, wx23,  wy23;
real[] t24, wl24, wu24, wx24,  wy24;
real[] t25, wl25, wu25, wx25,  wy25;
real[] t26, wl26, wu26, wx26,  wy26;
real[] t27, wl27, wu27, wx27,  wy27;
real[] t28, wl28, wu28, wx28,  wy28;
real[] t29, wl29, wu29, wx29,  wy29;

real[] t30, wl30, wu30, wx30,  wy30;
real[] t30, wl30, wu30, wx30,  wy30;
real[] t31, wl31, wu31, wx31,  wy31;
real[] t32, wl32, wu32, wx32,  wy32;
real[] t33, wl33, wu33, wx33,  wy33;
real[] t34, wl34, wu34, wx34,  wy34;
real[] t35, wl35, wu35, wx35,  wy35;
real[] t36, wl36, wu36, wx36,  wy36;
real[] t37, wl37, wu37, wx37,  wy37;
real[] t38, wl38, wu38, wx38,  wy38;
real[] t39, wl39, wu39, wx39,  wy39;

real[] t40, wl40, wu40, wx40,  wy40;
real[] t40, wl40, wu40, wx40,  wy40;
real[] t41, wl41, wu41, wx41,  wy41;
real[] t42, wl42, wu42, wx42,  wy42;
real[] t43, wl43, wu43, wx43,  wy43;
real[] t44, wl44, wu44, wx44,  wy44;
real[] t45, wl45, wu45, wx45,  wy45;
real[] t46, wl46, wu46, wx46,  wy46;
real[] t47, wl47, wu47, wx47,  wy47;
real[] t48, wl48, wu48, wx48,  wy48;
real[] t49, wl49, wu49, wx49,  wy49;

real[] t50, wl50, wu50, wx50,  wy50;
real[] t50, wl50, wu50, wx50,  wy50;
real[] t51, wl51, wu51, wx51,  wy51;
real[] t52, wl52, wu52, wx52,  wy52;
real[] t53, wl53, wu53, wx53,  wy53;
real[] t54, wl54, wu54, wx54,  wy54;
real[] t55, wl55, wu55, wx55,  wy55;
real[] t56, wl56, wu56, wx56,  wy56;
real[] t57, wl57, wu57, wx57,  wy57;
real[] t58, wl58, wu58, wx58,  wy58;
real[] t59, wl59, wu59, wx59,  wy59;

for (int j = 0; j < N; ++j)
  {
   if ((int) m[j] == 1)
      {
	t1.push (tt[j]);
 	wl1.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu1.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
        wx1.push (psi[j]-w[j]/2.);
	wy1.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 2)
      {
	t2.push (tt[j]);
  	wl2.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu2.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx2.push (psi[j]-w[j]/2.);
	wy2.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 3)
      {
	t3.push (tt[j]);
 	wl3.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu3.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx3.push (psi[j]-w[j]/2.);
	wy3.push (psi[j]+w[j]/2.);	
       }
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
   if ((int) m[j] == 17)
      {
	t17.push (tt[j]);
 	wl17.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu17.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx17.push (psi[j]-w[j]/2.);
	wy17.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 18)
      {
	t18.push (tt[j]);
     	wl18.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu18.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx18.push (psi[j]-w[j]/2.);
	wy18.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 19)
      {
	t19.push (tt[j]);
	wl19.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu19.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx19.push (psi[j]-w[j]/2.);
	wy19.push (psi[j]+w[j]/2.);	
      }

   if ((int) m[j] == 20)
      {
	t20.push (tt[j]);
	wl20.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu20.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx20.push (psi[j]-w[j]/2.);
	wy20.push (psi[j]+w[j]/2.);	
      }
 if ((int) m[j] == 21)
      {
	t21.push (tt[j]);
 	wl21.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu21.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
        wx21.push (psi[j]-w[j]/2.);
	wy21.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 22)
      {
	t22.push (tt[j]);
  	wl22.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu22.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx22.push (psi[j]-w[j]/2.);
	wy22.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 23)
      {
	t23.push (tt[j]);
 	wl23.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu23.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx23.push (psi[j]-w[j]/2.);
	wy23.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 24)
      {
	t24.push (tt[j]);
     	wl24.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu24.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx24.push (psi[j]-w[j]/2.);
	wy24.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 25)
      {
	t25.push (tt[j]);
	wl25.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu25.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx25.push (psi[j]-w[j]/2.);
	wy25.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 26)
      {
	t26.push (tt[j]);
	wl26.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu26.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx26.push (psi[j]-w[j]/2.);
	wy26.push (psi[j]+w[j]/2.);	
      }
   if ((int) m[j] == 27)
      {
	t27.push (tt[j]);
 	wl27.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu27.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx27.push (psi[j]-w[j]/2.);
	wy27.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 28)
      {
	t28.push (tt[j]);
     	wl28.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu28.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx28.push (psi[j]-w[j]/2.);
	wy28.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 29)
      {
	t29.push (tt[j]);
	wl29.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu29.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx29.push (psi[j]-w[j]/2.);
	wy29.push (psi[j]+w[j]/2.);	
      }

if ((int) m[j] == 30)
      {
	t30.push (tt[j]);
	wl30.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu30.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx30.push (psi[j]-w[j]/2.);
	wy30.push (psi[j]+w[j]/2.);	
      }
 if ((int) m[j] == 31)
      {
	t31.push (tt[j]);
 	wl31.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu31.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
        wx31.push (psi[j]-w[j]/2.);
	wy31.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 32)
      {
	t32.push (tt[j]);
  	wl32.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu32.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx32.push (psi[j]-w[j]/2.);
	wy32.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 33)
      {
	t33.push (tt[j]);
 	wl33.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu33.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx33.push (psi[j]-w[j]/2.);
	wy33.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 34)
      {
	t34.push (tt[j]);
     	wl34.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu34.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx34.push (psi[j]-w[j]/2.);
	wy34.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 35)
      {
	t35.push (tt[j]);
	wl35.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu35.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx35.push (psi[j]-w[j]/2.);
	wy35.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 36)
      {
	t36.push (tt[j]);
	wl36.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu36.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx36.push (psi[j]-w[j]/2.);
	wy36.push (psi[j]+w[j]/2.);	
      }
   if ((int) m[j] == 37)
      {
	t37.push (tt[j]);
 	wl37.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu37.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx37.push (psi[j]-w[j]/2.);
	wy37.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 38)
      {
	t38.push (tt[j]);
     	wl38.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu38.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx38.push (psi[j]-w[j]/2.);
	wy38.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 39)
      {
	t39.push (tt[j]);
	wl39.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu39.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx39.push (psi[j]-w[j]/2.);
	wy39.push (psi[j]+w[j]/2.);	
      }

if ((int) m[j] == 40)
      {
	t40.push (tt[j]);
	wl40.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu40.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx40.push (psi[j]-w[j]/2.);
	wy40.push (psi[j]+w[j]/2.);	
      }
 if ((int) m[j] == 41)
      {
	t41.push (tt[j]);
 	wl41.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu41.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
        wx41.push (psi[j]-w[j]/2.);
	wy41.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 42)
      {
	t42.push (tt[j]);
  	wl42.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu42.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx42.push (psi[j]-w[j]/2.);
	wy42.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 43)
      {
	t43.push (tt[j]);
 	wl43.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu43.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx43.push (psi[j]-w[j]/2.);
	wy43.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 44)
      {
	t44.push (tt[j]);
     	wl44.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu44.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx44.push (psi[j]-w[j]/2.);
	wy44.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 45)
      {
	t45.push (tt[j]);
	wl45.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu45.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx45.push (psi[j]-w[j]/2.);
	wy45.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 46)
      {
	t46.push (tt[j]);
	wl46.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu46.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx46.push (psi[j]-w[j]/2.);
	wy46.push (psi[j]+w[j]/2.);	
      }
   if ((int) m[j] == 47)
      {
	t47.push (tt[j]);
 	wl47.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu47.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx47.push (psi[j]-w[j]/2.);
	wy47.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 48)
      {
	t48.push (tt[j]);
     	wl48.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu48.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx48.push (psi[j]-w[j]/2.);
	wy48.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 49)
      {
	t49.push (tt[j]);
	wl49.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu49.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx49.push (psi[j]-w[j]/2.);
	wy49.push (psi[j]+w[j]/2.);	
      }

if ((int) m[j] == 50)
      {
	t50.push (tt[j]);
	wl50.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu50.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx50.push (psi[j]-w[j]/2.);
	wy50.push (psi[j]+w[j]/2.);	
      }
 if ((int) m[j] == 51)
      {
	t51.push (tt[j]);
 	wl51.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu51.push (psi[j]+cos(tt[j]*pi)*w[j]/2.); 
        wx51.push (psi[j]-w[j]/2.);
	wy51.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 52)
      {
	t52.push (tt[j]);
  	wl52.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu52.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx52.push (psi[j]-w[j]/2.);
	wy52.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 53)
      {
	t53.push (tt[j]);
 	wl53.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu53.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx53.push (psi[j]-w[j]/2.);
	wy53.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 54)
      {
	t54.push (tt[j]);
     	wl54.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu54.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx54.push (psi[j]-w[j]/2.);
	wy54.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 55)
      {
	t55.push (tt[j]);
	wl55.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu55.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx55.push (psi[j]-w[j]/2.);
	wy55.push (psi[j]+w[j]/2.);	
      }
    if ((int) m[j] == 56)
      {
	t56.push (tt[j]);
	wl56.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu56.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx56.push (psi[j]-w[j]/2.);
	wy56.push (psi[j]+w[j]/2.);	
      }
   if ((int) m[j] == 57)
      {
	t57.push (tt[j]);
 	wl57.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu57.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx57.push (psi[j]-w[j]/2.);
	wy57.push (psi[j]+w[j]/2.);	
       }
    if ((int) m[j] == 58)
      {
	t58.push (tt[j]);
     	wl58.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu58.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx58.push (psi[j]-w[j]/2.);
	wy58.push (psi[j]+w[j]/2.);	
 	}
    if ((int) m[j] == 59)
      {
	t59.push (tt[j]);
	wl59.push (psi[j]-cos(tt[j]*pi)*w[j]/2.);
	wu59.push (psi[j]+cos(tt[j]*pi)*w[j]/2.);
        wx59.push (psi[j]-w[j]/2.);
	wy59.push (psi[j]+w[j]/2.);	
      }
  }

pen s;

if (t1.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t1, wl1), s);
     draw (graph (t1, wu1), s);
     draw (graph (t1, wx1), s);
     draw (graph (t1, wy1), s);
   }
if (t2.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t2, wl2), s);
     draw (graph (t2, wu2), s);
     draw (graph (t2, wx2), s);
     draw (graph (t2, wy2), s);
    }	
if (t3.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t3, wl3), s);
     draw (graph (t3, wu3), s);
     draw (graph (t3, wx3), s);
     draw (graph (t3, wy3), s);
   }	
if (t4.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t4, wl4), s);
     draw (graph (t4, wu4), s);
     draw (graph (t4, wx4), s);
     draw (graph (t4, wy4), s);
   }
if (t5.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t5, wl5), s);
     draw (graph (t5, wu5), s);
     draw (graph (t5, wx5), s);
     draw (graph (t5, wy5), s);
  } 	
if (t6.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t6, wl6), s);
     draw (graph (t6, wu6), s);
     draw (graph (t6, wx6), s);
     draw (graph (t6, wy6), s);
   }	
if (t7.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t7, wl7), s);
     draw (graph (t7, wu7), s);
     draw (graph (t7, wx7), s);
     draw (graph (t7, wy7), s);
   }
if (t8.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t8, wl8), s);
     draw (graph (t8, wu8), s);
     draw (graph (t8, wx8), s);
     draw (graph (t8, wy8), s);
   }
if (t9.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t9, wl9), s);
     draw (graph (t9, wu9), s);
     draw (graph (t9, wx9), s);
     draw (graph (t9, wy9), s);
   }

if (t10.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t10, wl10), s);
     draw (graph (t10, wu10), s);
     draw (graph (t10, wx10), s);
     draw (graph (t10, wy10), s);
   }
if (t11.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t11, wl11), s);
     draw (graph (t11, wu11), s);
     draw (graph (t11, wx11), s);
     draw (graph (t11, wy11), s);
   }
if (t12.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t12, wl12), s);
     draw (graph (t12, wu12), s);
     draw (graph (t12, wx12), s);
     draw (graph (t12, wy12), s);
    }	
if (t13.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t13, wl13), s);
     draw (graph (t13, wu13), s);
     draw (graph (t13, wx13), s);
     draw (graph (t13, wy13), s);
   }	
if (t14.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t14, wl14), s);
     draw (graph (t14, wu14), s);
     draw (graph (t14, wx14), s);
     draw (graph (t14, wy14), s);
   }
if (t15.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t15, wl15), s);
     draw (graph (t15, wu15), s);
     draw (graph (t15, wx15), s);
     draw (graph (t15, wy15), s);
  } 	
if (t16.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t16, wl16), s);
     draw (graph (t16, wu16), s);
     draw (graph (t16, wx16), s);
     draw (graph (t16, wy16), s);
   }	
if (t17.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t17, wl17), s);
     draw (graph (t17, wu17), s);
     draw (graph (t17, wx17), s);
     draw (graph (t17, wy17), s);
   }
if (t18.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t18, wl18), s);
     draw (graph (t18, wu18), s);
     draw (graph (t18, wx18), s);
     draw (graph (t18, wy18), s);
   }
if (t19.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t19, wl19), s);
     draw (graph (t19, wu19), s);
     draw (graph (t19, wx19), s);
     draw (graph (t19, wy19), s);
   }

if (t20.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t20, wl20), s);
     draw (graph (t20, wu20), s);
     draw (graph (t20, wx20), s);
     draw (graph (t20, wy20), s);
   }
  if (t21.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t21, wl21), s);
     draw (graph (t21, wu21), s);
     draw (graph (t21, wx21), s);
     draw (graph (t21, wy21), s);
   }
if (t22.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t22, wl22), s);
     draw (graph (t22, wu22), s);
     draw (graph (t22, wx22), s);
     draw (graph (t22, wy22), s);
    }	
if (t23.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t23, wl23), s);
     draw (graph (t23, wu23), s);
     draw (graph (t23, wx23), s);
     draw (graph (t23, wy23), s);
   }	
if (t24.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t24, wl24), s);
     draw (graph (t24, wu24), s);
     draw (graph (t24, wx24), s);
     draw (graph (t24, wy24), s);
   }
if (t25.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t25, wl25), s);
     draw (graph (t25, wu25), s);
     draw (graph (t25, wx25), s);
     draw (graph (t25, wy25), s);
  } 	
if (t26.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t26, wl26), s);
     draw (graph (t26, wu26), s);
     draw (graph (t26, wx26), s);
     draw (graph (t26, wy26), s);
   }	
if (t27.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t27, wl27), s);
     draw (graph (t27, wu27), s);
     draw (graph (t27, wx27), s);
     draw (graph (t27, wy27), s);
   }
if (t28.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t28, wl28), s);
     draw (graph (t28, wu28), s);
     draw (graph (t28, wx28), s);
     draw (graph (t28, wy28), s);
   }
if (t29.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t29, wl29), s);
     draw (graph (t29, wu29), s);
     draw (graph (t29, wx29), s);
     draw (graph (t29, wy29), s);
   }

if (t30.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t30, wl30), s);
     draw (graph (t30, wu30), s);
     draw (graph (t30, wx30), s);
     draw (graph (t30, wy30), s);
   }
  if (t31.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t31, wl31), s);
     draw (graph (t31, wu31), s);
     draw (graph (t31, wx31), s);
     draw (graph (t31, wy31), s);
   }
if (t32.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t32, wl32), s);
     draw (graph (t32, wu32), s);
     draw (graph (t32, wx32), s);
     draw (graph (t32, wy32), s);
    }	
if (t33.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t33, wl33), s);
     draw (graph (t33, wu33), s);
     draw (graph (t33, wx33), s);
     draw (graph (t33, wy33), s);
   }	
if (t34.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t34, wl34), s);
     draw (graph (t34, wu34), s);
     draw (graph (t34, wx34), s);
     draw (graph (t34, wy34), s);
   }
if (t35.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t35, wl35), s);
     draw (graph (t35, wu35), s);
     draw (graph (t35, wx35), s);
     draw (graph (t35, wy35), s);
  } 	
if (t36.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t36, wl36), s);
     draw (graph (t36, wu36), s);
     draw (graph (t36, wx36), s);
     draw (graph (t36, wy36), s);
   }	
if (t37.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t37, wl37), s);
     draw (graph (t37, wu37), s);
     draw (graph (t37, wx37), s);
     draw (graph (t37, wy37), s);
   }
if (t38.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t38, wl38), s);
     draw (graph (t38, wu38), s);
     draw (graph (t38, wx38), s);
     draw (graph (t38, wy38), s);
   }
if (t39.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t39, wl39), s);
     draw (graph (t39, wu39), s);
     draw (graph (t39, wx39), s);
     draw (graph (t39, wy39), s);
   }

if (t40.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t40, wl40), s);
     draw (graph (t40, wu40), s);
     draw (graph (t40, wx40), s);
     draw (graph (t40, wy40), s);
   }
  if (t41.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t41, wl41), s);
     draw (graph (t41, wu41), s);
     draw (graph (t41, wx41), s);
     draw (graph (t41, wy41), s);
   }
if (t42.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t42, wl42), s);
     draw (graph (t42, wu42), s);
     draw (graph (t42, wx42), s);
     draw (graph (t42, wy42), s);
    }	
if (t43.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t43, wl43), s);
     draw (graph (t43, wu43), s);
     draw (graph (t43, wx43), s);
     draw (graph (t43, wy43), s);
   }	
if (t44.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t44, wl44), s);
     draw (graph (t44, wu44), s);
     draw (graph (t44, wx44), s);
     draw (graph (t44, wy44), s);
   }
if (t45.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t45, wl45), s);
     draw (graph (t45, wu45), s);
     draw (graph (t45, wx45), s);
     draw (graph (t45, wy45), s);
  } 	
if (t46.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t46, wl46), s);
     draw (graph (t46, wu46), s);
     draw (graph (t46, wx46), s);
     draw (graph (t46, wy46), s);
   }	
if (t47.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t47, wl47), s);
     draw (graph (t47, wu47), s);
     draw (graph (t47, wx47), s);
     draw (graph (t47, wy47), s);
   }
if (t48.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t48, wl48), s);
     draw (graph (t48, wu48), s);
     draw (graph (t48, wx48), s);
     draw (graph (t48, wy48), s);
   }
if (t49.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t49, wl49), s);
     draw (graph (t49, wu49), s);
     draw (graph (t49, wx49), s);
     draw (graph (t49, wy49), s);
   }

if (t50.length > 0)
   {
     s = black + solid + 0.5;
     draw (graph (t50, wl50), s);
     draw (graph (t50, wu50), s);
     draw (graph (t50, wx50), s);
     draw (graph (t50, wy50), s);
   }
  if (t51.length > 0)
   {
     s = red + solid + 0.5;
     draw (graph (t51, wl51), s);
     draw (graph (t51, wu51), s);
     draw (graph (t51, wx51), s);
     draw (graph (t51, wy51), s);
   }
if (t52.length > 0)
   {
     s = yellow + solid + 0.5;
     draw (graph (t52, wl52), s);
     draw (graph (t52, wu52), s);
     draw (graph (t52, wx52), s);
     draw (graph (t52, wy52), s);
    }	
if (t53.length > 0)
   {
     s = green + solid + 0.5;
     draw (graph (t53, wl53), s);
     draw (graph (t53, wu53), s);
     draw (graph (t53, wx53), s);
     draw (graph (t53, wy53), s);
   }	
if (t54.length > 0)
   {
     s = blue + solid + 0.5;
     draw (graph (t54, wl54), s);
     draw (graph (t54, wu54), s);
     draw (graph (t54, wx54), s);
     draw (graph (t54, wy54), s);
   }
if (t55.length > 0)
   {
     s = orange + solid + 0.5;
     draw (graph (t55, wl55), s);
     draw (graph (t55, wu55), s);
     draw (graph (t55, wx55), s);
     draw (graph (t55, wy55), s);
  } 	
if (t56.length > 0)
   {
     s = purple + solid + 0.5;
     draw (graph (t56, wl56), s);
     draw (graph (t56, wu56), s);
     draw (graph (t56, wx56), s);
     draw (graph (t56, wy56), s);
   }	
if (t57.length > 0)
   {
     s = magenta + solid + 0.5;
     draw (graph (t57, wl57), s);
     draw (graph (t57, wu57), s);
     draw (graph (t57, wx57), s);
     draw (graph (t57, wy57), s);
   }
if (t58.length > 0)
   {
     s = cyan + solid + 0.5;
     draw (graph (t58, wl58), s);
     draw (graph (t58, wu58), s);
     draw (graph (t58, wx58), s);
     draw (graph (t58, wy58), s);
   }
if (t59.length > 0)
   {
     s = pink + solid + 0.5;
     draw (graph (t59, wl59), s);
     draw (graph (t59, wu59), s);
     draw (graph (t59, wx59), s);
     draw (graph (t59, wy59), s);
   }
   
s = dotted + black + 1.5;
ylimits (0.85, 1., Crop);
yequals (psilim, s);
yequals (psiped, s);

pen qq = fontsize (25.);
defaultpen (qq);
xaxis ("$t ({\rm ms})$", BottomTop, LeftTicks);
yaxis ("${\mit\Psi}_N$", LeftRight, RightTicks);
