# RESCALE

## Description
    Rescales equilibrium gFile to change q_95 by modifying toroidal plasma current 
    while keeping B_toroidal the same

## Requirements

   - asymptote (for plots) (https://asymptote.sourceforge.io)	
   
## Contents	

### /Documentation
- Rescale.tex: Latex description of program algorithm

### /Source 
- Makefile: GNU makefile 
- *.f90: Fortran_90 source files
 
### /Inputs
- gFile: Initial gFile
- q_95: Contains target q_95
	
### /Outputs
- gFile: Rescaled gFile
 	
### /Plots 
- *.asy: Asymptote scripts to plot rescaled equilibrium

## gFile Format

    read (100, '(a48, 3i4)') string, i3, NRBOX, NZBOX
    read (100, '(5e16.9)') RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero
    read (100, '(5e16.9)') RAXIS, ZAXIS, PSIAXIS, PSIBOUND, B0
    read (100, '(5e16.9)') CURRENT, zero, zero, zero, zero
    read (100, '(5e16.9)') zero, zero, zero, zero, zero
 
    read (100, '(5e16.9)') (T (i), i = 1, NRBOX)
    read (100, '(5e16.9)') (P (i), i = 1, NRBOX)
    read (100, '(5e16.9)') (TTp (i), i = 1, NRBOX)
    read (100, '(5e16.9)') (Pp (i), i = 1, NRBOX)
  
    read (100, '(5e16.9)') ((PSI (i, j), i = 1, NRBOX), j = 1, NZBOX)

    read (100, '(5e16.9)') (Q (i), i = 1, NRBOX)

    read (100, '(2i5)') NBOUND, NLIM

    read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
    read (100, '(5e16.9)') (RLIM (i), ZLIM (i), i = 1, NLIM)
  
  *Anything after this is ignored*
  
