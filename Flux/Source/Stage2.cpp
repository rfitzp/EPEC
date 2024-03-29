// Stage2.cpp

// PROGRAM ORGANIATION:
//
// void Flux:: Stage2                      ()
// void Flux:: WriteStage2Netcdfcpp        ()
// void Flux:: WriteStage2Netcdfc          ()
// void Flux:: Stage2SetSimpsonWeights     ()
// void Flux:: Stage2ReadData              ()
// void Flux:: Stage2CalcQ                 ()
// void Flux:: Stage2FindRational          ()
// void Flux:: Stage2CalcGGJ               ()
// void Flux:: Stage2CalcStraightAngle     ()
// void Flux:: Stage2CalcNeoclassicalAngle ()
// void Flux:: Stage2CalcNeoclassicalPara  ()
// void Flux:: Stage2CalcMatrices          ()
// void Flux:: Stage2CalcTearing           ()

#include "Flux.h"

// ####################################################
// Function to input Stage1 data and output Stage2 data
// ####################################################
void Flux::Stage2 ()
{
  // ..............................
  // Set weights for Simpson's rule
  // ..............................
  Stage2SetSimpsonWeights ();
  
  // ................................
  // Read data for Stage2 calculation
  // ................................
  Stage2ReadData ();
  fflush (stdout);
  
  // ..........................
  // Calculate Stage2 q profile
  // ..........................
  Stage2CalcQ ();
  fflush (stdout);

  // ......................
  // Find rational surfaces
  // ......................
  Stage2FindRational ();
  fflush (stdout);

  // .......................................
  // Calculate GGJ data at rational surfaces
  // .......................................
  Stage2CalcGGJ ();
  fflush (stdout);
  
  // ..................................................
  // Calculate straight angle data on rational surfaces
  // ..................................................
  Stage2CalcStraightAngle ();
  fflush (stdout);

  // ......................................................
  // Calculate neoclassical angle data on rational surfaces
  // ......................................................
  Stage2CalcNeoclassicalAngle ();
  fflush (stdout);
 
  // ......................................................
  // Calculate neoclassical parameters at rational surfaces
  // ......................................................
  Stage2CalcNeoclassicalPara ();
  fflush (stdout);
  
  // ....................................
  // Calculate tearing stability matrices
  // ....................................
  Stage2CalcMatrices ();
  fflush (stdout);

  // ...............................................
  // Calculate cylindrical tearing stability indices
  // ...............................................
  Stage2CalcTearing ();
  fflush (stdout);

  // ..................
  // Output NETCDF file
  // ..................
  WriteStage2Netcdfc ();
 
  // ............
  // Output fFile
  // ............
  FILE* file = OpenFilew ((char*) "Outputs/fFile");
  fprintf (file, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %d %d %d %16.9e %16.9e %16.9e %16.6e\n",
	   R0, B0, ra * R0, q95, r95 /ra, qrat, rrat /ra, QP[0], QP[NPSI-1], NPSI, NTOR, nres, PSILIM, PSIPED, Pped, PSIRAT);
  for (int j = 0; j < NPSI; j++)
    fprintf (file, "%16.9e %16.9e %16.9e\n",
	     1. - P[j], rP[j] /ra, - ra * Interpolate (NPSI, rP, P, rP[j], 1));
  for (int i = 0; i < nres; i++)
    fprintf (file, "%d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	     mres[i], rres[i]/ra, sres[i], gres[i], gmres[i], Ktres[i], Kares[i], fcres[i], ajj[i], PsiNres[i], dPsidr[i], Khres[i], A1res[i], A2res[i], q_hat[i], C1res[i], C2res[i], E[i]+F[i]+H[i]*H[i], Poem1[i], Poem2[i], Poem3[i], Delta_w[i], Sigma_w[i]);
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      fprintf (file, "%d %d %16.9e %16.9e\n", mres[i], mres[j],
	       GSL_REAL (gsl_matrix_complex_get (FF, i, j)), GSL_IMAG (gsl_matrix_complex_get (FF, i, j)));
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      fprintf (file, "%d %d %16.9e %16.9e\n", mres[i], mres[j],
	       GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
  fclose (file);
  
  // ........
  // Clean up
  // ........
  delete[] RPTS;  delete[] ZPTS;
  delete[] RBPTS; delete[] ZBPTS;
  delete[] RLPTS; delete[] ZLPTS;

  delete[] PSIN; delete[] G;   delete[] Pr;
  delete[] GGp;  delete[] Prp; delete[] Q;
  delete[] Jphi;

  delete[] s; delete[] Rs;

  delete[] P;    delete[] RP;  delete[] rP; delete[] GP;
  delete[] QGP;  delete[] QP;  delete[] PP; delete[] GPP;
  delete[] PPP;  delete[] S;   delete[] QX; delete[] RP1;
  delete[] Bt;   delete[] Bt1; delete[] Bp; delete[] Bp1;
  delete[] J0;   

  delete[] PsiN; delete[] QPN; delete[] QPPN, delete[] A1, delete[] A2;
  
  delete[] mres;    delete[] qres;  delete[] rres;  delete[] sres;
  delete[] gres;    delete[] Rres;  delete[] gmres; delete[] fcres;
  delete[] PsiNres; delete[] Rres1; delete[] A1res; delete[] A2res;
  delete[] JP;      delete[] JPr;   delete[] JPrr;  delete[] alpres;
  delete[] betres;  delete[] gamres;

  delete[] th; 
  gsl_matrix_free (Rst); gsl_matrix_free (Zst);

  delete[] Th;
  gsl_matrix_free (Rnc); gsl_matrix_free (Znc); 

  gsl_matrix_free (Bnc);    gsl_matrix_free (Cnc); gsl_matrix_free (theta);
  gsl_matrix_free (factor);
 
  delete[] I1;    delete[] I2;    delete[] I3;    delete[] I7;    delete[] I8;  
  delete[] Ktres; delete[] Kares; delete[] ajj;   delete[] dPsidr;
  delete[] Khres; delete[] q_hat; delete[] C1res; delete[] C2res;

  delete[] J1; delete[] J2;  delete[] J3; delete[] J4; delete[] J5;
  delete[] J6; delete[] E;   delete[] F;  delete[] H;

  gsl_matrix_free (I4); gsl_matrix_free (I5); gsl_matrix_free (I6);
  
  gsl_matrix_complex_free (FF); gsl_matrix_complex_free (EE);

  delete[] A;        delete[] B;       delete[] Delta_nw; delete[] Sigma_nw;
  delete[] Poem1;    delete[] Poem2;   delete[] Poem3;    delete[] Delta_pw;
  delete[] Sigma_pw; delete[] Delta_w; delete[] Sigma_w;
}

// ####################################
// Function to write Stage2 NETCDF file
// ####################################
void Flux::WriteStage2Netcdfc ()
{
  // Open file
  int err = 0, dataFile;
  err = nc_create ("Outputs/Stage2.nc", NC_CLOBBER, &dataFile);
  
  if (err != 0)
    {
      printf ("FLUX::WriteStage2Netcdfc: Error opening Outputs/Stage2.nc\n");
      exit (1);
    }
  
  // Basic parameters
  int para_d, para;
  double parameters[12] = {R0, B0, RLEFT, ZLOW, RRIGHT, ZHIGH, Raxis, Zaxis, q95, PSILIM, PSIRAT, PSIPED};
  err += nc_def_dim (dataFile, "index", 12, &para_d);
  err += nc_def_var (dataFile, "Parameters", NC_DOUBLE, 1, &para_d, &para);
  nc_put_att_text   (dataFile, para, "Contents", strlen ("R0 B0 RLEFT ZLOW RRIGH ZHIGH RAXIS ZAXIS q95 PSILIM PSIRAT PSIPED"),
		     "R0 B0 RLEFT ZLOW RRIGH ZHIGH RAXIS ZAXIS q95 PSILIM PSIRAT PSIPED");
  
  // Plasma boundary
  int bound_d, bound_r, bound_z;
  err += nc_def_dim (dataFile, "i_bound", NBPTS, &bound_d);
  err += nc_def_var (dataFile, "RBPTS", NC_DOUBLE, 1, &bound_d, &bound_r);
  err += nc_def_var (dataFile, "ZBPTS", NC_DOUBLE, 1, &bound_d, &bound_z);
  
  // R
  int R_d, R;
  err += nc_def_dim (dataFile, "i", NRPTS, &R_d);
  err += nc_def_var (dataFile, "R", NC_DOUBLE, 1, &R_d, &R);
 
  // Z
  int Z_d, Z;
  err += nc_def_dim (dataFile, "j", NZPTS, &Z_d);
  err += nc_def_var (dataFile, "Z", NC_DOUBLE, 1, &Z_d, &Z);
 
  // Psi
  double* DATA = new double[NRPTS*NZPTS];
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      DATA[j + i*NZPTS] = PSIARRAY (i, j);

  int psi_d[2], psi;
  psi_d[0] = R_d;
  psi_d[1] = Z_d;
  err += nc_def_var (dataFile, "PSI", NC_DOUBLE, 2, psi_d, &psi);
  
  // Psi_R
  double* DATA1 = new double[NRPTS*NZPTS];
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      DATA1[j + i*NZPTS] = GetPsiR (RPTS[i], ZPTS[j]);
  
  int psi_r;
  err += nc_def_var (dataFile, "PSI_R", NC_DOUBLE, 2, psi_d, &psi_r);
  
  // Psi_Z
  double* DATA2 = new double[NRPTS*NZPTS];
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      DATA2[j + i*NZPTS] = GetPsiZ (RPTS[i], ZPTS[j]);

  int psi_z;
  err += nc_def_var (dataFile, "PSI_Z", NC_DOUBLE, 2, psi_d, &psi_z);
 
  // Psi_N
  int Q_d, PN;
  err += nc_def_dim (dataFile, "N_PSI", NPSI, &Q_d);
  err += nc_def_var (dataFile, "PsiN", NC_DOUBLE, 1, &Q_d, &PN);

  // q
  int Q;
  err += nc_def_var (dataFile, "q", NC_DOUBLE, 1, &Q_d, &Q);

  // q_psi
  int Qp;
  err += nc_def_var (dataFile, "q_psi", NC_DOUBLE, 1, &Q_d, &Qp);
 
  // g
  int G;
  err += nc_def_var (dataFile, "g", NC_DOUBLE, 1, &Q_d, &G);

  // P
  int Px;
  err += nc_def_var (dataFile, "P", NC_DOUBLE, 1, &Q_d, &Px);
 
  // g_psi
  int Gp;
  err += nc_def_var (dataFile, "g_psi", NC_DOUBLE, 1, &Q_d, &Gp);
  
  // P_psi
  int Pp;
  err += nc_def_var (dataFile, "P_psi", NC_DOUBLE, 1, &Q_d, &Pp);

  // r/a
  double* rra = new double[NPSI];
  for (int j = 0; j < NPSI; j++)
    rra [j] = rP[j] /ra;
  int r_y;
  err += nc_def_var (dataFile, "r", NC_DOUBLE, 1, &Q_d, &r_y);

  // J
  int J_y;
  err += nc_def_var (dataFile, "J", NC_DOUBLE, 1, &Q_d, &J_y);

  // dJdr
  int dJdr_y;
  err += nc_def_var (dataFile, "dJdr", NC_DOUBLE, 1, &Q_d, &dJdr_y);    

  // d2Jdr2
  int d2Jdr2_y;
  err += nc_def_var (dataFile, "d2Jdr2", NC_DOUBLE, 1, &Q_d, &d2Jdr2_y);    
    
  // m_res
  int S_d, mr;
  err += nc_def_dim (dataFile, "i_res", nres, &S_d);
  err += nc_def_var (dataFile, "m_pol", NC_INT, 1, &S_d, &mr);

  // PsiN_res
  int PNr;
  err += nc_def_var (dataFile, "PsiN_res", NC_DOUBLE, 1, &S_d, &PNr);
  
  // q_res
  int qr;
  err += nc_def_var (dataFile, "q_res", NC_DOUBLE, 1, &S_d, &qr);
  
  // s_res
  int sr;
  err += nc_def_var (dataFile, "s_res", NC_DOUBLE, 1, &S_d, &sr);

  // g_res
  int gr;
  err += nc_def_var (dataFile, "g_res", NC_DOUBLE, 1, &S_d, &gr);
 
  // gamma_res
  int gmr;
  err += nc_def_var (dataFile, "gamma_res", NC_DOUBLE, 1, &S_d, &gmr);
  
  // f_c
  int fcr;
  err += nc_def_var (dataFile, "fc_res", NC_DOUBLE, 1, &S_d, &fcr);

  // a_jj
  int ajr; 
  err += nc_def_var (dataFile, "ajj_res", NC_DOUBLE, 1, &S_d, &ajr);
 
  // Q
  int qhr;
  err += nc_def_var (dataFile, "Q_res", NC_DOUBLE, 1, &S_d, &qhr);
 
  // A1
  int a1r;
  err +=  nc_def_var (dataFile, "A1_res", NC_DOUBLE, 1, &S_d, &a1r);

  // A2
  int a2r;
  err += nc_def_var (dataFile, "A2_res", NC_DOUBLE, 1, &S_d, &a2r);

  // D_R
  double* D_R = new double[nres];
  for (int i = 0; i < nres; i++)
    D_R[i] = E[i] + F[i] + H[i]*H[i];
  
  int dr;
  err += nc_def_var (dataFile, "DR_res", NC_DOUBLE, 1, &S_d, &dr);

  // Rst
  double* DATA3 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA3[j + i*NTHETA] = gsl_matrix_get (Rst, i, j);

  int rst_d[2], T_d, rst;
  err += nc_def_dim (dataFile, "k", NTHETA, &T_d);
  rst_d[0] = S_d;
  rst_d[1] = T_d;
  err += nc_def_var (dataFile, "R_st", NC_DOUBLE, 2, rst_d, &rst);

  // Zst
  double* DATA4 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA4[j + i*NTHETA] = gsl_matrix_get (Zst, i, j);

  int zst;
  err += nc_def_var (dataFile, "Z_st", NC_DOUBLE, 2, rst_d, &zst);

  // Rnc
  double* DATA5 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA5[j + i*NTHETA] = gsl_matrix_get (Rnc, i, j);

  int rnc;
  err += nc_def_var (dataFile, "R_nc", NC_DOUBLE, 2, rst_d, &rnc);

  // Znc
  double* DATA6 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA6[j + i*NTHETA] = gsl_matrix_get (Znc, i, j);

  int znc;
  err += nc_def_var (dataFile, "Z_nc", NC_DOUBLE, 2, rst_d, &znc);

  // theta
  double* DATA7 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA7[j + i*NTHETA] = gsl_matrix_get (theta, i, j);

  int the;
  err += nc_def_var (dataFile, "theta", NC_DOUBLE, 2, rst_d, &the);

  // B_nc
  double* DATA10 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA10[j + i*NTHETA] = gsl_matrix_get (Bnc, i, j);

  int bnc;
  err += nc_def_var (dataFile, "B_nc", NC_DOUBLE, 2, rst_d, &bnc);

  // B_nc
  double* DATA11 = new double[nres*NTHETA];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA11[j + i*NTHETA] = gsl_matrix_get (Cnc, i, j);

  int cnc;
  err += nc_def_var (dataFile, "C_nc", NC_DOUBLE, 2, rst_d, &cnc);

  // E_real
  double* DATA8 = new double[nres*nres];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      DATA8[j + i*nres] = GSL_REAL (gsl_matrix_complex_get (EE, i, j));

  int est_d[2], ereal;
  est_d[0] = S_d;
  est_d[1] = S_d;
  err += nc_def_var (dataFile, "E_real", NC_DOUBLE, 2, est_d, &ereal);

  // E_imag
  double* DATA9 = new double[nres*nres];
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      DATA9[j + i*nres] = GSL_IMAG (gsl_matrix_complex_get (EE, i, j));

  int eimag;
  err += nc_def_var (dataFile, "E_imag", NC_DOUBLE, 2, est_d, &eimag);

  // Delta_w
  int dw;
  err += nc_def_var (dataFile, "Delta_w", NC_DOUBLE, 1, &S_d, &dw);

  // Sigma_w
  int sw;
  err += nc_def_var (dataFile, "Sigma_w", NC_DOUBLE, 1, &S_d, &sw);

  // Poem1
  int p1;
  err += nc_def_var (dataFile, "Poem1", NC_DOUBLE, 1, &S_d, &p1);

  // Poem2
  int p2;
  err += nc_def_var (dataFile, "Poem2", NC_DOUBLE, 1, &S_d, &p2);

  // Poem3
  int p3;
  err += nc_def_var (dataFile, "Poem3", NC_DOUBLE, 1, &S_d, &p3);

  // A
  int aa;
  err += nc_def_var (dataFile, "A", NC_DOUBLE, 1, &S_d, &aa);

  // B
  int bb;
  err += nc_def_var (dataFile, "B", NC_DOUBLE, 1, &S_d, &bb);
  
  err += nc_enddef (dataFile);

  if (err != 0)
    {
      printf ("FLUX::WriteStage2Netcdfc: Error defining variables in Outputs/Stage2.nc\n");
      exit (1);
    }

  // Write data
  err += nc_put_var_double (dataFile, para,     parameters);
  err += nc_put_var_double (dataFile, bound_r,  RBPTS);
  err += nc_put_var_double (dataFile, bound_z,  ZBPTS);
  err += nc_put_var_double (dataFile, R,        RPTS);
  err += nc_put_var_double (dataFile, Z,        ZPTS);
  err += nc_put_var_double (dataFile, psi,      DATA);
  err += nc_put_var_double (dataFile, psi_r,    DATA1);
  err += nc_put_var_double (dataFile, psi_z,    DATA2);
  err += nc_put_var_double (dataFile, PN,       P);
  err += nc_put_var_double (dataFile, Q,        QP);
  err += nc_put_var_double (dataFile, Qp,       QPN);
  err += nc_put_var_double (dataFile, G,        GP);
  err += nc_put_var_double (dataFile, Gp,       GPP);
  err += nc_put_var_double (dataFile, Px,       PP);
  err += nc_put_var_double (dataFile, Pp,       PPP);
  err += nc_put_var_double (dataFile, r_y,      rra);
  err += nc_put_var_double (dataFile, J_y,      JP);
  err += nc_put_var_double (dataFile, dJdr_y,   JPr);
  err += nc_put_var_double (dataFile, d2Jdr2_y, JPrr);
  err += nc_put_var_int    (dataFile, mr,       mres);
  err += nc_put_var_double (dataFile, PNr,      PsiNres);
  err += nc_put_var_double (dataFile, qr,       qres);
  err += nc_put_var_double (dataFile, sr,       sres);
  err += nc_put_var_double (dataFile, gr,       gres);
  err += nc_put_var_double (dataFile, gmr,      gmres);
  err += nc_put_var_double (dataFile, fcr,      fcres);
  err += nc_put_var_double (dataFile, ajr,      ajj);
  err += nc_put_var_double (dataFile, qhr,      q_hat);
  err += nc_put_var_double (dataFile, a1r,      A1res);
  err += nc_put_var_double (dataFile, a2r,      A2res);
  err += nc_put_var_double (dataFile, dr,       D_R);
  err += nc_put_var_double (dataFile, rst,      DATA3);
  err += nc_put_var_double (dataFile, zst,      DATA4);
  err += nc_put_var_double (dataFile, rnc,      DATA5);
  err += nc_put_var_double (dataFile, znc,      DATA6);
  err += nc_put_var_double (dataFile, the,      DATA7);
  err += nc_put_var_double (dataFile, bnc,      DATA10);
  err += nc_put_var_double (dataFile, cnc,      DATA11);
  err += nc_put_var_double (dataFile, ereal,    DATA8);
  err += nc_put_var_double (dataFile, eimag,    DATA9);
  err += nc_put_var_double (dataFile, dw,       Delta_w);
  err += nc_put_var_double (dataFile, sw,       Sigma_w);
  err += nc_put_var_double (dataFile, p1,       Poem1);
  err += nc_put_var_double (dataFile, p2,       Poem2);
  err += nc_put_var_double (dataFile, p3,       Poem3);
  err += nc_put_var_double (dataFile, aa,       A);
  err += nc_put_var_double (dataFile, bb,       B);
  
  if (err != 0)
    {
      printf ("FLUX::WriteStage2Netcdfc: Error writing Outputs/Stage2.nc\n");
      exit (1);
    }
  
  // Close file
  err += nc_close (dataFile);

  if (err != 0)
    {
      printf ("FLUX::WriteStage2Netcdfc: Error closing Outputs/Stage2.nc\n");
      exit (1);
    }
   
  delete[] DATA;   delete[] DATA1, delete[] DATA2, delete[] D_R;
  delete[] DATA3;  delete[] DATA4; delete[] DATA5; delete[] DATA6;
  delete[] DATA7;  delete[] DATA8; delete[] DATA9; delete[] DATA10;
  delete[] DATA11; delete[] rra;
 }

// #############################################
// Function to assign weights for Simpson's rule
// #############################################
void Flux::Stage2SetSimpsonWeights ()
{
  double h  = 2.*M_PI   /double (NTHETA-1);  // Step length for angular integrals
         hh = 0.9999999 /double (NTHETA-1);  // Step length for fraction of circulating particles integrals

  Weight1D.resize (NTHETA);
  weight1D.resize (NTHETA);
  Weight2D.resize (NTHETA, NTHETA);
  for (int j = 0; j < NTHETA; j++)
    {
    	Weight1D (j) = 0.;
	weight1D (j) = 0.;

	for (int k = 0; k < NTHETA; k++)
	  Weight2D (j, k) = 0.;
    }
  for (int j = 0; j < NTHETA-2; j += 2)
    {
      Weight1D (j)   +=      h /3.;
      Weight1D (j+1) += 4. * h /3.;
      Weight1D (j+2) +=      h /3.;

      weight1D (j)   +=      hh /3.;
      weight1D (j+1) += 4. * hh /3.;
      weight1D (j+2) +=      hh /3.;
      
      for (int k = 0; k < NTHETA-2; k += 2)
	{
	  Weight2D (j,   k)   +=       h*h /9.;
	  Weight2D (j+1, k)   +=  4. * h*h /9.;
	  Weight2D (j+2, k)   +=       h*h /9.;
	  Weight2D (j,   k+1) +=  4. * h*h /9.;
	  Weight2D (j+1, k+1) += 16. * h*h /9.;
	  Weight2D (j+2, k+1) +=  4. * h*h /9.;
	  Weight2D (j,   k+2) +=       h*h /9.;
	  Weight2D (j+1, k+2) +=  4. * h*h /9.;
	  Weight2D (j+2, k+2) +=       h*h /9.;
	}
    }
}

// #############################################
// Function to read data for Stage2 calculations
// #############################################
void Flux::Stage2ReadData ()
{
  printf ("Reading data from gFile:\n");
  fflush (stdout);
  
  // ..............
  // Read R0 and B0
  // ..............
  FILE* file = OpenFiler ((char*) "Outputs/Stage1/R0B0.txt");
  if (fscanf (file, "%lf %lf", &R0, &B0) != 2)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/R0B0.txt\n");
      exit (1);
    }
  fclose (file);

  // ......................
  // Read bounding box data
  // ......................
  file = OpenFiler ((char*) "Outputs/Stage1/Box.txt");
  if (fscanf (file, "%lf %lf %lf %lf", &RLEFT, &ZLOW, &RRIGHT, &ZHIGH) != 4)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Box.txt\n");
      exit (1);
    }
  fclose (file);

  // .......................
  // Read magnetic axis data
  // .......................
  file = OpenFiler ((char*) "Outputs/Stage1/Axis.txt");
  if (fscanf (file, "%lf %lf", &Raxis, &Zaxis) != 2)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Axis.txt\n");
      exit (1);
    }
  fclose (file);

  // ................
  // Read array sizes
  // ................
  file = OpenFiler ((char*) "Outputs/Stage1/Points.txt");
  if (fscanf (file, "%d %d %d %d", &NRPTS, &NZPTS, &NBPTS, &NLPTS) != 4)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Points.txt\n");
      exit (1);
    }
  fclose (file);

  // ..............
  // Read R, Z grid
  // ..............
  RPTS = new double[NRPTS];  // R array
  ZPTS = new double[NZPTS];  // Z array

  file = OpenFiler ((char*) "Outputs/Stage1/R.txt");
  for (int i = 0; i < NRPTS; i++)
    if (fscanf (file, "%lf", &RPTS[i]) != 1)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/R.txt\n");
	exit (1);
      }
  fclose (file);
  file = OpenFiler ((char*) "Outputs/Stage1/Z.txt");
  for (int j = 0; j < NZPTS; j++)
    if (fscanf (file, "%lf", &ZPTS[j]) != 1)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Z.txt\n");
	exit (1);
      }
  fclose (file);

  dR  = RPTS[1] - RPTS[0];
  dZ  = ZPTS[1] - ZPTS[0];
  dR2 = dR*dR;
  dZ2 = dZ*dZ;
  dR3 = dR*dR2;
  dZ3 = dZ*dZ2;

  // ..............................
  // Read boundary and limiter data
  // ..............................
  RBPTS = new double[NBPTS];  // R coordinates of boundary
  ZBPTS = new double[NBPTS];  // Z coordinates of boundary

  file = OpenFiler ((char*) "Outputs/Stage1/Boundary.txt");
  for (int i = 0; i < NBPTS; i++)
    if (fscanf (file, "%lf %lf", &RBPTS[i], &ZBPTS[i]) != 2)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Boundary.txt\n");
	exit (1);
      }
  fclose (file);

  RLPTS = new double[NLPTS];  // R coordinates of limiter
  ZLPTS = new double[NLPTS];  // Z coordinates of limiter

  file = OpenFiler ((char*) "Outputs/Stage1/Limiter.txt");
  for (int i = 0; i < NLPTS; i++)
    if (fscanf (file, "%lf %lf", &RLPTS[i], &ZLPTS[i]) != 2)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Limiter.txt\n");
	exit (1);
      }
  fclose (file);
  printf ("R0      = %11.4e  B0          = %11.4e\n", R0,  B0);

  // ........
  // Read Psi
  // ........
  PSIARRAY.resize (NRPTS, NZPTS);  // Psi (R, Z)

  double val; int ival;
  file = OpenFiler ((char*) "Outputs/Stage1/PsiSequential.txt");
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {	
	if (fscanf (file, "%d %d %lf", &ival, &ival, &val) != 3)
	  {
	    printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/PsiSequential.txt\n");
	    exit (1);
	  }
	PSIARRAY (i, j) = val;
      }
  fclose (file);

  // ....................................................
  // Modify PSI such that PsiBoundary = 0 and PsiAxis = 1
  // ....................................................
  double PsiAxis     = InterpolatePsi (Raxis,    Zaxis,    0);
  double PsiBoundary = InterpolatePsi (RBPTS[0], ZBPTS[0], 0);

  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {
	double val = PSIARRAY (i, j) - PsiBoundary;
	PSIARRAY (i, j) = val;
      }

  PsiAxis     = InterpolatePsi (Raxis,    Zaxis,    0);
  PsiBoundary = InterpolatePsi (RBPTS[0], ZBPTS[0], 0);
  Psic        = PsiAxis / (R0*R0*B0);

  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {
	double val = PSIARRAY (i, j) /PsiAxis;
	PSIARRAY (i, j) = val;
      }

  printf ("PsiAxis = %11.4e  PsiBoundary = %11.4e  PsiAxis /(R0*R0*B0) = %11.4e\n", PsiAxis, PsiBoundary, Psic);

  // .............................
  // Read equilibrium profile data
  // .............................
  PSIN = new double[NRPTS]; // PSI_N array 
  G    = new double[NRPTS]; // g
  Pr   = new double[NRPTS]; // p
  GGp  = new double[NRPTS]; // g dg/dpsi
  Prp  = new double[NRPTS]; // dp/dpsi
  Q    = new double[NRPTS]; // q
  Jphi = new double[NRPTS]; // J_phi
  double dum;

  file = OpenFiler ((char*) "Outputs/Stage1/Profiles.txt");
  for (int i = 0; i < NRPTS; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf %lf", &PSIN[i], &G[i], &Pr[i], &GGp[i], &Prp[i], &Q[i], &Jphi[i]) != 7)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Stage1/Profiles.txt\n");
	exit (1);
      }
  fclose (file);

  // ....................................................
  // Find closest equilibrium grid-point to magnetic axis
  // ....................................................
  double rmin = 1.e6;
  for (int i = 0; i < NRPTS; i++)
    if (fabs (RPTS[i] - Raxis) < rmin)
      {
	rmin = fabs (RPTS[i] - Raxis);
	ic   = i;
      }
  rmin = 1.e6;
  for (int j = 0; j < NZPTS; j++)
    if (fabs (ZPTS[j] - Zaxis) < rmin)
      {
	rmin = fabs (ZPTS[j] - Zaxis);
	jc   = j;
      }

  // ..........................................................................
  // Find closest equilibrium grid-point to inner magnetic boundary on midplane
  // ..........................................................................
  ia = 0;
  for (int i = 0; i <= ic; i++)
    {
      if (PSIARRAY (i, jc) < 0.)
	ia++;
    }
  Rbound = (RPTS[ia-1]*PSIARRAY (ia, jc) - RPTS[ia]*PSIARRAY (ia-1, jc))
    /(PSIARRAY (ia, jc) - PSIARRAY (ia-1, jc));

  L = ic - ia + 2; // Number of points in Psi (R, Zaxis) array

  // ..........................................................................
  // Find closest equilibrium grid-point to outer magnetic boundary on midplane
  // ..........................................................................
  int ia1 = ic;
  for (int i = ic; i < NRPTS; i++)
    {
      if (PSIARRAY (i, jc) > 0.)
	ia1++;
    }
  Rbound1 = (RPTS[ia1+1]*PSIARRAY (ia1, jc) - RPTS[ia1]*PSIARRAY (ia1+1, jc))
    /(PSIARRAY (ia1, jc) - PSIARRAY (ia1+1, jc));
}

// ######################################
// Function to calculate Stage2 q profile
// ######################################
void Flux::Stage2CalcQ ()
{
  // ............................
  // Setup Psi (R, Zaxis) profile
  // ............................
  s = new double[L]; // Array of s = sqrt [1 - Psi (R, Zaxis)] values

  for (int l = 0; l < L-2; l++)
    s[L-2-l] = sqrt (1. - PSIARRAY (l+ia, jc));
  s[0]   = 0.;
  s[L-1] = 1.;

  // ...........................
  // Setup R (Z = Zaxis) profile
  // ...........................
  Rs = new double[L]; // Array of R(s) values

  for (int l = 0; l < L-2; l++)
    Rs[L-2-l] = RPTS[l+ia];
  Rs[0]   = Raxis;
  Rs[L-1] = Rbound;

  // ......................
  // Set up Stage2 Psi grid
  // ......................
  P   = new double[NPSI];  // 1 - Psi array
  RP  = new double[NPSI];  // R(Psi) on inboard midplane
  RP1 = new double[NPSI];  // R(Psi) on outboard midplane
  Bt  = new double[NPSI];  // B_toroidal(Psi) on inboard midplane
  Bt1 = new double[NPSI];  // B_toroidal(Psi) on outboard midplane
  Bp  = new double[NPSI];  // B_poloidal(Psi) on inboard midplane
  Bp1 = new double[NPSI];  // B_poloidal(Psi) on outboard midplane
  rP  = new double[NPSI];  // r(Psi)
  GP  = new double[NPSI];  // g(Psi)
  QGP = new double[NPSI];  // q(Psi)/g(Psi) 
  QP  = new double[NPSI];  // q(Psi)
  PP  = new double[NPSI];  // P(Psi)
  GPP = new double[NPSI];  // dg/dPsi
  PPP = new double[NPSI];  // dP/dPsi
  S   = new double[NPSI];  // sqrt (1 - Psi)
  QX  = new double[NPSI];  // q(Psi) from gFile
  J0  = new double[NPSI];  // GGJ integral

  PsiN = new double[NPSI]; // PsiN array
  QPN  = new double[NPSI]; // dQ/dPsiN array
  QPPN = new double[NPSI]; // d^2Q/dPsiN^2 array
  A1   = new double[NPSI]; // QP/QPN/fabs(Psic) array
  A2   = new double[NPSI]; // QPPN/QPN/3 array
  JP   = new double[NPSI]; // J_phi
  JPr  = new double[NPSI]; // dJ_phi/dr
  JPrr = new double[NPSI]; // d^2J_phi/dr^2
  
  for (int j = 0; j < NPSI; j++)
    {
      double s = double (j) /double (NPSI-1);

      PsiN[j] = PSILIM * (1. - pow (1. - s, PACK));
      P   [j] = 1. - PsiN[j];
      S   [j] = sqrt (1. - P[j]);
    }

  // .......................................
  // Calculate Stage2 g(Psi), P(Psi) profile
  // .......................................
  for (int j = 0; j < NPSI; j++)
    {
      GP [j] = Interpolate (NRPTS, PSIN, G,    PsiN[j], 0);
      PP [j] = Interpolate (NRPTS, PSIN, Pr,   PsiN[j], 0);
      GPP[j] = Interpolate (NRPTS, PSIN, GGp,  PsiN[j], 0) /GP[j];
      PPP[j] = Interpolate (NRPTS, PSIN, Prp,  PsiN[j], 0);
      QX [j] = Interpolate (NRPTS, PSIN, Q,    PsiN[j], 0);
      JP [j] = Interpolate (NRPTS, PSIN, Jphi, PsiN[j], 0);
    }

  // ...............................
  // Calculate Stage2 R(Psi) profile
  // ...............................
  for (int j = 0; j < NPSI; j++)
    RP[j] = Interpolate (L, s, Rs, S[j], 0);
  
  // ......................................
  // Calculate Stage2 q(Psi)/g(Psi) profile
  // ......................................
  printf ("Calculating q(Psi)/g(Psi) profile:\n");
  fflush (stdout);
  CalcQGP ();
  
  RP1[0]      = Raxis;
  RP1[NPSI-1] = Rbound1;

  QGP[0] = Q[0] /G[0];
  QP [0] = Q[0];

  // Calculate and smooth QPN profile (quartic interpolation)
  double* QPSmooth = new double[NPSI];

  for (int i = 0; i < NPSI; i++)
    QPSmooth[i] = QP[i];

  for (int i = 0; i < NSMOOTH; i++)
    Smoothing (NPSI, QPSmooth);

  for (int j = 0; j < NPSI; j++)
    {
      QPN  [j] = InterpolateQ (NPSI, PsiN, QPSmooth, PsiN[j], 1, 4);
      QPPN [j] = InterpolateQ (NPSI, PsiN, QPSmooth, PsiN[j], 2, 4);
    }
   
  for (int i = 0; i < NSMOOTH; i++)
    {
      Smoothing (NPSI, QPN);
      Smoothing (NPSI, QPPN);
    }

  delete[] QPSmooth;
  
  // Calculate A1 profile
  for (int j = 0; j < NPSI; j++)
    {
      A1[j] = QP[j] /QPN[j] /fabs(Psic);
      A2[j] = QPPN[j] /QPN[j] /3.;
    }

  for (int i = 0; i < NSMOOTH; i++)
    {
      Smoothing (NPSI, A1);
      Smoothing (NPSI, A2);
    }

  // ...........................................
  // Calculate B_tor, B_pol profiles on midplane
  // ...........................................
  for (int j = 0; j < NPSI; j++)
    {
      Bt[j]  = GP[j] /RP[j];
      Bt1[j] = GP[j] /RP1[j];

      double PsiR = GetPsiR (RP[j], Zaxis);
      double PsiZ = GetPsiZ (RP[j], Zaxis);

      Bp[j] = fabs(Psic) * sqrt (PsiR*PsiR + PsiZ*PsiZ) /RP[j];

      PsiR = GetPsiR (RP1[j], Zaxis);
      PsiZ = GetPsiZ (RP1[j], Zaxis);

      Bp1[j] = fabs(Psic) * sqrt (PsiR*PsiR + PsiZ*PsiZ) /RP1[j];
    }

  // .............................
  // Calculate Stage2 r(P) profile
  // .............................
  printf ("Calculating r(Psi) profile:\n");
  fflush (stdout);
  CalcrP ();
  rP[0] = 0.;

  // ...............................................
  // Calculate and smooth Jp, Jpr, and Jprr profiles
  // ...............................................
  for (int i = 0; i < NSMOOTH; i++)
    Smoothing (NPSI, JP);
  
  for (int j = 0; j < NPSI; j++)
    JPr [j] = Interpolate (NPSI, rP, JP, rP[j], 1);
 
  for (int i = 0; i < NSMOOTH; i++)
    Smoothing (NPSI, JPr);
  
  for (int j = 0; j < NPSI; j++)
    JPrr[j] = Interpolate (NPSI, rP, JPr, rP[j], 1);
    
  for (int i = 0; i < NSMOOTH; i++)
    Smoothing (NPSI, JPrr);
 
  // .................................
  // Find q95, r95, qrat, rrat, qa, ra
  // .................................
  double g95;
  double s95  = sqrt (0.95);
  double sped = sqrt (PSIPED);
  double srat = sqrt (PSIRAT);
  q95    = Interpolate (NPSI, S, QP, s95,  0);
  r95    = Interpolate (NPSI, S, rP, s95,  0);
  g95    = Interpolate (NPSI, S, GP, s95,  0);
  Pped   = Interpolate (NPSI, S, PP, sped, 0)/ PP[0];
  qrat   = Interpolate (NPSI, S, QP, srat,  0);
  rrat   = Interpolate (NPSI, S, rP, srat,  0);
  qa     = QP[NPSI-1];
  ra     = rP[NPSI-1];
  printf ("q95 = %9.2e  r95/ra = %9.2e  qrat = %9.2e  rrat/ra = %9.2e  qa = %9.2e  ra = %9.2e  a = %9.2e (m)  Pped/P(0) = %9.2e\n",
	  q95, r95 /ra, qrat, rrat/ra, qa, ra, ra*R0, Pped);
}

// ###################################
// Function to find rational surfaces.
// Assumes monotonic q-profile.
// ###################################
void Flux::Stage2FindRational ()
{
  // .....................................
  // Determine number of rational surfaces
  // .....................................
  int mmin = int (double (NTOR) * QP[0]) + 1;
  int mmax = int (double (NTOR) * qrat);
  if (mmin < MMIN)
    mmin = MMIN;
  if (mmax > MMAX)
    mmax = MMAX;
  nres = mmax - mmin + 1;

  if (nres <= 0)
    {
      printf ("FLUX::Stage2FindRational: Error no rational surfaces in plasma\n");
      exit (1);
    }

  // .....................................
  // Calculate rational surface quantities
  // .....................................
  mres    = new int   [nres];
  qres    = new double[nres];
  rres    = new double[nres];
  sres    = new double[nres];
  gres    = new double[nres];
  Rres    = new double[nres];
  Rres1   = new double[nres];
  gmres   = new double[nres];
  fcres   = new double[nres];
  PsiNres = new double[nres];
  A1res   = new double[nres];
  A2res   = new double[nres];
  alpres  = new double[nres];
  betres  = new double[nres];
  gamres  = new double[nres];

  for (int i = 0; i < nres; i++)
    {
      mres[i] = mmin + i;

      qres[i] = double (mres[i]) /double (NTOR);
      rres[i] = Interpolate (NPSI, QP, rP, qres[i], 0);
      sres[i] = rres[i] /Interpolate (NPSI, QP, rP, qres[i], 1) /qres[i];

      // Check that PsiN < PSILIM
      PsiNres[i] = 1. - Interpolate (NPSI, rP, P, rres[i], 0);
      if (PsiNres[i] > PSILIM)
	{
	  nres = i;
	  break; 
	}

      // Correct rres values
      for (int ii = 0; ii < 4; ii++)
	{
	  double qqq = Interpolate (NPSI, rP, QP, rres[i], 0);
	  double qqp = Interpolate (NPSI, rP, QP, rres[i], 1);
	  double qpp = Interpolate (NPSI, rP, QP, rres[i], 2);
	  rres[i] += (- qqp + sqrt (qqp*qqp - 2.*qpp * (qqq - qres[i]))) /qpp;
	}
      
      gmres  [i] = Interpolate (NPSI, rP, QP, rres[i], 0) - qres[i];
      sres   [i] = rres[i] /Interpolate (NPSI, QP, rP, qres[i], 1) /qres[i];
      gres   [i] = Interpolate (NPSI, rP, GP,   rres[i], 0);
      Rres   [i] = Interpolate (NPSI, rP, RP,   rres[i], 0);
      Rres1  [i] = Interpolate (NPSI, rP, RP1,  rres[i], 0);
      PsiNres[i] = Interpolate (NPSI, rP, PsiN, rres[i], 0);
      A1res  [i] = Interpolate (NPSI, rP, A1,   rres[i], 0);
      A2res  [i] = Interpolate (NPSI, rP, A2,   rres[i], 0);

      alpres [i] = - Interpolate (NPSI, rP, JPr, rres[i], 0) * qres[i] /sres[i];
      betres [i] = alpres[i] * (Interpolate (NPSI, rP, JPrr, rres[i], 0) /Interpolate (NPSI, rP, JPr, rres[i], 0)
				+ (sres[i] - 1.) /rres[i]
				- Interpolate (NPSI, rP, QP, rres[i], 2) /2. /Interpolate (NPSI, rP, QP, rres[i], 1));
      gamres [i] = - Interpolate (NPSI, rP, JPrr, rres[i], 0) * qres[i] /sres[i];
    }

  printf ("Rational surface data:\n");
  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d  PsiNs = %9.2e  rs/ra = %9.2e  ss = %9.2e  residual = %9.2e  R = %9.2e  A1 = %9.2e  residual = %9.2e  A2 = %9.2e  alp = %9.2e  bet = %9.2e  gam = %9.2e\n",
	    mres[i], PsiNres[i], rres[i] /ra, sres[i], gmres[i], Rres1[i], A1res[i],
	    1. - (gres[i]/sres[i]/qres[i]) *rres[i]*rres[i] /Psic/Psic /A1res[i], A2res[i], alpres[i], betres[i], gamres[i]);
   
  // .....................................
  // Confirm q values at rational surfaces
  // .....................................
  printf ("Confirm rational surface q-values:\n");
  CheckQP ();
}

// ##########################################################
// Calculate Glasser-Greene-Johnson data at rational surfaces
// ##########################################################
void Flux::Stage2CalcGGJ ()
{
  J1 = new double[nres];
  J2 = new double[nres];
  J3 = new double[nres];
  J4 = new double[nres];
  J5 = new double[nres];
  J6 = new double[nres];
  E  = new double[nres];
  F  = new double[nres];
  H  = new double[nres];

  printf ("Calculating Glasser-Greene-Johnson data at rational surfaces:\n");
  fflush (stdout);
  CalcGGJ ();

  // Calculate E, F, H
  for (int i = 0; i < nres; i++)
    {
      double qpsi = Interpolate (NPSI, rP, QP, rres[i], 1);
      double ppsi = Interpolate (NPSI, rP, PP, rres[i], 1);
      double J1p  = Interpolate (NPSI, rP, J0, rres[i], 1);
      
      E[i] = - (ppsi /qpsi/qpsi)      * (J1p - gres[i] * qpsi * J1[i] /J2[i]) * J5[i];
      F[i] =   (ppsi*ppsi /qpsi/qpsi) * (gres[i]*gres[i] * (J5[i] * J6[i] - J4[i] * J4[i]) + J3[i] * J5[i]);
      H[i] =   (ppsi /qpsi)           * (J4[i] - J1[i] * J5[i] /J2[i]) * gres[i];
    }

  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d J1 = %9.2e J2 = %9.2e J3 = %9.2e J4 = %9.2e J5 = %9.2e J6 = %9.2e E = %9.2e F = %10.4e H = %9.2e DI = %9.2e DR_ln = %9.2e DR_nl = %9.2e\n",
	    mres[i], J1[i], J2[i], J3[i], J4[i], J5[i], J6[i], E[i], F[i], H[i], E[i]+F[i]+H[i]-0.25, E[i]+F[i]+H[i]*H[i], E[i]+F[i]);
}

// ##################################################
// Calculate straight angle data at rational surfaces
// ##################################################
void Flux::Stage2CalcStraightAngle ()
{
  // ..................
  // Set up theta array
  // ..................
  th = new double[NTHETA]; 
  for (int k = 0; k < NTHETA; k++)
    {
      double t = double (k) /double (NTHETA-1);
      th[k]    = 2.* M_PI * t;
    }
  Rst = gsl_matrix_alloc (nres, NTHETA);
  Zst = gsl_matrix_alloc (nres, NTHETA);

  // .............................
  // Calculate straight angle data
  // .............................
  printf ("Calculating straight angle data at rational surfaces:\n");
  fflush (stdout);
  CalcStraightAngle ();
}

// ######################################################
// Calculate neoclassical angle data at rational surfaces
// ######################################################
void Flux::Stage2CalcNeoclassicalAngle ()
{
  // ...........................................
  // Calculate gamma values at rational surfaces
  // ...........................................
  printf ("Calculating gamma values at rational surfaces:\n");
  fflush (stdout);
  CalcGamma ();

  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d  rs/ra = %9.2e  g = %9.2e  gamma = %9.2e  gamma*q/g = %9.2e\n",
	    mres[i], rres[i] /ra, gres[i], gmres[i], gmres[i]*qres[i] /gres[i]);

  // ..................
  // Set up Theta array
  // ..................
  Th = new double[NTHETA];
  for (int k = 0; k < NTHETA; k++)
    {
      double t = double (k) /double (NTHETA-1);
      Th[k]    = 2.* M_PI * t;
    }
  Rnc    = gsl_matrix_alloc (nres, NTHETA);
  Znc    = gsl_matrix_alloc (nres, NTHETA);
  theta  = gsl_matrix_alloc (nres, NTHETA);
  factor = gsl_matrix_alloc (nres, NTHETA);

  // .................................
  // Calculate neoclassical angle data
  // .................................
  printf ("Calculating neoclassical angle data at rational surfaces:\n");
  CalcNeoclassicalAngle ();

  // Calculate transformation factor between theta and Theta
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < NTHETA; j++)
      {
	double R    = gsl_matrix_get (Rnc, i,  j);
	double Z    = gsl_matrix_get (Znc, i,  j);

	double PsiR = GetPsiR (R, Z);
	double PsiZ = GetPsiZ (R, Z);
	double Grad = PsiR*PsiR + PsiZ*PsiZ;
	double RB   = sqrt (gres[i]*gres[i] + Psic*Psic * Grad);

	double fac  = gres[i] /qres[i] /gmres[i] /RB /R;

	gsl_matrix_set (factor, i, j, fac);
      }

  // ..................................
  // Calculate |B| on rational surfaces
  // ..................................
  Bnc = gsl_matrix_alloc (nres, NTHETA);
  for (int j = 0; j < nres; j++)
    for (int k = 0; k < NTHETA; k++)
      {
	double gval = gres[j];
	double Rval = gsl_matrix_get (Rnc, j, k);
	double Zval = gsl_matrix_get (Znc, j, k);
	double PsiR = GetPsiR (Rval, Zval);
	double PsiZ = GetPsiZ (Rval, Zval);
	double Grad = sqrt (PsiR*PsiR + PsiZ*PsiZ);
	double Bval = sqrt (gval*gval + Psic*Psic * Grad*Grad) /Rval;

	gsl_matrix_set (Bnc, j, k, Bval);
      }

  // .........................................
  // Calculate d|B|dTheta on rational surfaces
  // .........................................
  Cnc = gsl_matrix_alloc (nres, NTHETA);
  double* Carray = new double[NTHETA];
  for (int j = 0; j < nres; j++)
    {
      for (int k = 0; k < NTHETA; k++)
	{
	  double gval  = gres [j];
	  double gmval = gmres[j];
	  double Rval  = gsl_matrix_get (Rnc, j, k);
	  double Zval  = gsl_matrix_get (Znc, j, k);
	  double PsiR  = GetPsiR  (Rval, Zval);
	  double PsiZ  = GetPsiZ  (Rval, Zval);
	  double PsiRR = GetPsiRR (Rval, Zval);
	  double PsiRZ = GetPsiRZ (Rval, Zval);
	  double PsiZZ = GetPsiZZ (Rval, Zval);
	  double Grad  = sqrt (PsiR*PsiR + PsiZ*PsiZ);
	  double Bval  = sqrt (gval*gval + Psic*Psic * Grad*Grad) /Rval;
	  double dRdT  = - fabs (Psic) * PsiZ /gmval /Bval /Rval;
	  double dZdT  = + fabs (Psic) * PsiR /gmval /Bval /Rval;
	  
	  double val   = - Bval*dRdT/Rval + (Psic*Psic/Bval/Rval/Rval) * ((PsiR*PsiRR + PsiZ*PsiRZ)*dRdT + (PsiR*PsiRZ + PsiZ*PsiZZ)*dZdT);

	  Carray[k] = val;
	}

      // Smooth data
      for (int i = 0; i < 10; i++)
	Smoothing (NTHETA, Carray);
	  
      for (int k = 0; k < NTHETA; k++)
	gsl_matrix_set (Cnc, j, k, Carray[k]);
    }
  delete[] Carray;
}

// ######################################################
// Calculate neoclassical parameters at rational surfaces
// ######################################################
void Flux::Stage2CalcNeoclassicalPara ()
{
  // ................
  // Allocate memory
  // ...............
  I1     = new double[nres];
  I2     = new double[nres];
  I3     = new double[nres];
  I4     = gsl_matrix_alloc (nres, NNC);
  I5     = gsl_matrix_alloc (nres, NNC);
  I6     = gsl_matrix_alloc (nres, NTHETA);
  I7     = new double[nres];
  I8     = new double[nres];
  Ktres  = new double[nres];
  Kares  = new double[nres];
  Khres  = new double[nres];
  ajj    = new double[nres];
  dPsidr = new double[nres];
  q_hat  = new double[nres];
  C1res  = new double[nres];
  C2res  = new double[nres];
  
  // ......................................................
  // Calculate neoclassical parameters at rational surfaces 
  // ......................................................
  printf ("Calculating neoclassical parameters at rational surfaces:\n");
  fflush (stdout);
  double sum;
  for (int i = 0; i < nres; i++)
    {
      // Calculate I1
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) /gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I1[i] = sum;

      // Calculate I2
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) * gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I2[i] = sum;

      // Calculate I3
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) * gsl_matrix_get (Cnc, i, j) * gsl_matrix_get (Cnc, i, j)
	  /gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I3[i] = sum;

      // Calculate I4
      for (int k = 0; k < NNC; k++)
	{
	  double kk = double (k+1);

	  sum = 0.;
	  for (int j = 0; j < NTHETA; j++)
	    sum += Weight1D (j) * cos (kk * Th[j]) /gsl_matrix_get (Bnc, i, j);
	  sum *= sqrt (2.*kk) /2./M_PI;

	  gsl_matrix_set (I4, i, k, sum);
	}

      // Calculate I5
      for (int k = 0; k < NNC; k++)
	{
	  double kk = double (k+1);

	  sum = 0.;
	  for (int j = 0; j < NTHETA; j++)
	    sum += Weight1D (j) * cos (kk * Th[j]) /gsl_matrix_get (Bnc, i, j)
	      /gsl_matrix_get (Bnc, i, j) /2.;
	  sum *= sqrt (2.*kk) /2./M_PI;

	  gsl_matrix_set (I5, i, k, sum);
	}

      // Calcuate Kt, Ka, and Kh
      sum = 0.;
      for (int k = 0; k < NNC; k++)
	sum += gsl_matrix_get (I4, i, k) * gsl_matrix_get (I5, i, k);
      Ktres[i] = I1[i]*I1[i] * I3[i] /I2[i]/I2[i] /sum;
      Kares[i] = (8./3./M_PI) * (I2[i] /I3[i]) * Ktres[i]*Ktres[i];
      Khres[i] = gres[i]*gres[i] * I1[i] /Rres1[i]/Rres1[i] /I2[i];

      // Calculate Bmax
      double Bmax = -1.;
      for (int k = 0; k < NTHETA; k++)
	if (gsl_matrix_get (Bnc, i, k) > Bmax)
	  Bmax = gsl_matrix_get (Bnc, i, k);
      
      // Calculate I6
      for (int k = 0; k < NTHETA; k++)
	{
	  double lambda = double (k) * hh;

	  sum = 0.;
	  for (int kk = 0; kk < NTHETA; kk++)
	    sum += Weight1D (kk) * sqrt (1. - lambda * gsl_matrix_get (Bnc, i, kk) /Bmax)
	      /gsl_matrix_get (Bnc, i, kk);
	  sum /= 2.*M_PI;
	  
	  gsl_matrix_set (I6, i, k, sum);
	}

      // Calculate I7
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) * gsl_matrix_get (Rnc, i, j) * gsl_matrix_get (Rnc, i, j)
	  /gsl_matrix_get (Bnc, i, j);
      sum /= 2.*M_PI;
      
      I7[i] = sum;

      // Calculate I8
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	sum += Weight1D (j) /gsl_matrix_get (Bnc, i, j) /gsl_matrix_get (Bnc, i, j) /gsl_matrix_get (Bnc, i, j)
	  /gsl_matrix_get (Rnc, i, j) /gsl_matrix_get (Rnc, i, j);
      sum /= 2.*M_PI;
      
      I8[i] = sum;

      // Calculate fraction of circulating particles
      sum = 0.;
      for (int k = 0; k < NTHETA; k++)
	{
	  double lambda = double (k) * hh;
	  sum += weight1D (k) * lambda /gsl_matrix_get (I6, i, k);
	}
      sum *= 0.75 * I2[i] /Bmax/Bmax;

      fcres[i] = sum;

      // Calculate ajj
      sum = 0.;
      for (int j = 0; j < NTHETA; j++)
	{
	  double Rval = gsl_matrix_get (Rst, i, j);
	  double Zval = gsl_matrix_get (Zst, i, j);
 	  double PsiR = GetPsiR (Rval, Zval);
	  double PsiZ = GetPsiZ (Rval, Zval);
	  sum += Weight1D (j) /(PsiR*PsiR + PsiZ*PsiZ);
	}
      double fac = rres[i] * gres[i] /qres[i] /Psic;
      sum *= fac*fac /(2.*M_PI);
     
      ajj[i] = sum;

      // Calculate dPsidr 
      dPsidr[i] = rres[i] * gres[i] /qres[i] /fabs(Psic);

      // Calculate q_hat 
      q_hat[i] = (qres[i] /rres[i]) * sqrt ((1. - I1[i]*I1[i] * gres[i] /I7[i] /gmres[i] /qres[i]) /2./ajj[i]);
      if (isnan(q_hat[i]))
	{
	  printf ("PHASE::STAGE2 - Error q_hat[%2d] is NaN\n", i);
	  exit (1);
	}

      // Calculate C1
      C1res[i] = gmres[i] * qres[i] /I1[i] /gres[i];

      // Calculate C2
      C2res[i] = gres[i]*gres[i] * I8[i] /I1[i];
    }

  // ..............................
  // Output neoclassical parameters
  // ..............................
  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d I1 = %9.2e I2 = %9.2e I3 = %9.2e I4 = (%9.2e; %9.2e) I5 = (%9.2e; %9.2e) Kt = %9.2e Ka = %9.2e Kh = %9.2e fc = %9.2e\n",
	    mres[i], I1[i], I2[i], I3[i],
	    gsl_matrix_get (I4, i, 0), gsl_matrix_get (I4, i, NNC-1),
	    gsl_matrix_get (I5, i, 0), gsl_matrix_get (I5, i, NNC-1),
	    Ktres[i], Kares[i], Khres[i], fcres[i]);
  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d ajj = %9.2e dPsiNdr = %9.2e q_hat = %9.2e C1 = %9.2e C2 = %9.2e\n",
	    mres[i], ajj[i], dPsidr[i], q_hat[i], C1res[i], C2res[i]);
}

// ####################################
// Calculate tearing stability matrices
// ####################################
void Flux::Stage2CalcMatrices ()
{
  printf ("Calculating tearing stability matrices:\n");
  fflush (stdout);
  
  // ..................
  // Calculate F matrix
  // ..................
  FF = gsl_matrix_complex_alloc (nres, nres);
  EE = gsl_matrix_complex_alloc (nres, nres);
 
  if (NEOANG)
    {
      for (int i = 0; i < nres; i++)
	{
	  for (int j = 0; j < nres; j++)
	    {
	      double sumc = 0.;
	      for (int k = 0; k < NTHETA; k++)
		for (int kk = 0; kk < NTHETA; kk++)
		  sumc += Weight2D (k, kk) * GreenPlasmaCosNeo (i, k, j, kk);
	      sumc /= 2.*M_PI * 2.*M_PI;
	      
	      double sums = 0.;
	      for (int k = 0; k < NTHETA; k++)
		for (int kk = 0; kk < NTHETA; kk++)
		  sums += Weight2D (k, kk) * GreenPlasmaSinNeo (i, k, j, kk);
	      sums /= 2.*M_PI * 2.*M_PI;
	      
	      gsl_matrix_complex_set (FF, i, j, gsl_complex_rect (sumc, sums));
	    }
	  printf ("mpol = %3d\n", mres[i]);
	  fflush (stdout);
	}
    }
  else
    {
      for (int i = 0; i < nres; i++)
	{
	  for (int j = 0; j < nres; j++)
	    {
	      double sumc = 0.;
	      for (int k = 0; k < NTHETA; k++)
		for (int kk = 0; kk < NTHETA; kk++)
		  sumc += Weight2D (k, kk) * GreenPlasmaCos (i, k, j, kk);
	      sumc /= 2.*M_PI * 2.*M_PI;
	      
	      double sums = 0.;
	      for (int k = 0; k < NTHETA; k++)
		for (int kk = 0; kk < NTHETA; kk++)
		  sums += Weight2D (k, kk) * GreenPlasmaSin (i, k, j, kk);
	      sums /= 2.*M_PI * 2.*M_PI;
	      
	      gsl_matrix_complex_set (FF, i, j, gsl_complex_rect (sumc, sums));
	    }
	  printf ("mpol = %3d\n", mres[i]);
	  fflush (stdout);
	}
    }
      
  printf ("F-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%9.2e,%9.2e) ", GSL_REAL (gsl_matrix_complex_get (FF, i, j)), GSL_IMAG (gsl_matrix_complex_get (FF, i, j)));
      printf ("\n");
    }

  // ------------------
  // Calculate E-matrix
  // ------------------
  gsl_permutation*    px  = gsl_permutation_alloc (nres);
  int                 sss = 0;
  gsl_matrix_complex* FFF = gsl_matrix_complex_alloc (nres, nres);

  for (int i = 0; i < nres; i++)
    for (int j = 0; j < nres; j++)
      gsl_matrix_complex_set (FFF, i, j, gsl_matrix_complex_get (FF, i, j));

  gsl_linalg_complex_LU_decomp (FFF, px, &sss);
  gsl_linalg_complex_LU_invert (FFF, px, EE);

  gsl_permutation_free    (px);
  gsl_matrix_complex_free (FFF);

  printf ("E-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%8.1e,%8.1e) ", GSL_REAL (gsl_matrix_complex_get (EE, i, j)), GSL_IMAG (gsl_matrix_complex_get (EE, i, j)));
      printf ("\n");
    }
}

// ###############################################
// Calculate cylindrical tearing stability indices
// ###############################################
void Flux::Stage2CalcTearing ()
{
  printf ("Calculating cylindrical tearing data:\n");
  fflush (stdout);

  A        = new double[nres];
  B        = new double[nres];
  Delta_nw = new double[nres];
  Sigma_nw = new double[nres];
  Delta_pw = new double[nres];
  Sigma_pw = new double[nres];
  Delta_w  = new double[nres];
  Sigma_w  = new double[nres];
  Poem1    = new double[nres];
  Poem2    = new double[nres];
  Poem3    = new double[nres];

  for (int i = 0; i < nres; i++)
    {
      CalcTearingSolutionNoWall      (i);
      CalcTearingSolutionPerfectWall (i);
      Delta_w[i] = - Sigma_w[i]*Sigma_w[i] /(Delta_nw[i] - Delta_pw[i]);

      A[i]     = alpres[i] * rres[i];
      B[i]     = gamres[i] * rres[i] * rres[i];
      Poem1[i] = 0.41 * A[i]*A[i];
      Poem2[i] = 0.41 * (A[i]*A[i] * (log(1./rres[i]) + 4.85) - 0.5 * A[i] * Sigma_nw[i] - B[i] - 0.44 * A[i]);
      Poem3[i] = 0.80 * A[i]*A[i] - 0.27 * B[i] - 0.09 * A[i];
 
      printf ("mpol = %3d rs/a = %9.2e Deltanw = %9.2e Sigmanw = %9.2e Deltapw = %9.2e Sigmapw  = %9.2e Deltaw = %9.2e Sigmaw = %9.2e A = %9.2e B = %9.2e POEM1 = %9.2e POEM2 = %9.2e POEM3 = %9.2e\n",
	      mres[i], rres[i]/ra, Delta_nw[i], Sigma_nw[i], Delta_pw[i], Sigma_pw[i], Delta_w[i], Sigma_w[i], A[i], B[i], Poem1[i], Poem2[i], Poem3[i]);
      fflush (stdout);
    }
}
